/// Monte Carlo importance sampling engine.
///
/// Two simulation functions:
///
/// 1. `simulate_pedigree_probability_batch` — estimates P(hv), the average
///    probability of the observed haplotype configuration in the pedigree.
///    Traverses BFS from root; for each unknown individual, samples a mutated
///    haplotype from its parent using biased mutation; collects the importance
///    weight factor u/w and the product of "fixed" edge probabilities.
///
/// 2. `simulate_matching_haplotypes_batch` — estimates P(m(Hu)=x | hv), the
///    probability that x unknown individuals match the suspect.
///    Picks one "fixed" individual (weighted by picking_probabilities), assigns
///    the suspect haplotype to it, then simulates all other unknowns.
use crate::{
    Bias, BiasDirection, HaplotypeClass, Haplotype, MarkerSet, Pedigree, SimulationParameters,
};
use crate::simulation::bias::compute_biases_for_individual;
use crate::simulation::mutation::mutate_haplotype;
use crate::simulation::probability::get_edge_probability;
use crate::graph::{bfs_layers, descendants, most_recent_common_ancestor, shortest_path_length};
use crate::Result;
use rand::prelude::*;
use rust_decimal::Decimal;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// BatchResult
// ---------------------------------------------------------------------------

/// Running statistics for one Monte Carlo batch (many iterations).
#[derive(Debug, Clone)]
pub struct BatchResult {
    /// Sum of (probability * importance_weight) — numerator
    pub weighted_sum: Decimal,
    /// Sum of importance_weights — denominator
    pub weight_sum: Decimal,
    /// Number of iterations completed
    pub iterations: u64,
    /// Per-match-count weighted accumulators: match_count → weighted_sum
    pub match_accumulators: HashMap<u32, Decimal>,
    /// Per-individual weighted accumulators: individual_id → weighted_sum
    pub per_individual: HashMap<String, Decimal>,
}

impl BatchResult {
    pub fn new() -> Self {
        Self {
            weighted_sum: Decimal::ZERO,
            weight_sum: Decimal::ZERO,
            iterations: 0,
            match_accumulators: HashMap::new(),
            per_individual: HashMap::new(),
        }
    }

    /// Running importance-weighted mean estimate.
    pub fn running_mean(&self) -> Option<Decimal> {
        if self.weight_sum == Decimal::ZERO {
            return None;
        }
        Some(self.weighted_sum / self.weight_sum)
    }

    /// Add one iteration's result.
    pub fn accumulate(&mut self, probability: Decimal, importance_weight: Decimal) {
        self.weight_sum += importance_weight;
        self.weighted_sum += probability * importance_weight;
        self.iterations += 1;
    }
}

impl Default for BatchResult {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Bias computation (mirrors Python's Pedigree.get_biases)
// ---------------------------------------------------------------------------

/// Compute biases for an unknown individual in the context of a specific
/// haplotype configuration.
///
/// Matches Python's Pedigree.get_biases:
/// - Finds the parent of the individual from `child_of` map.
/// - Finds known descendants (individuals with known/suspect haplotype class).
/// - For each marker, checks if ALL known descendants agree on mutation direction
///   from parent to descendant; if so, creates a Bias.
/// - Uses distance_to_mrca for the default bias target mass.
fn get_biases_for_individual(
    pedigree: &Pedigree,
    individual_id: &str,
    marker_set: &MarkerSet,
    current_haplotypes: &HashMap<String, Haplotype>,
    bias_value: Option<f64>,
) -> Vec<Bias> {
    // bias_value <= 0 means disable bias entirely
    if let Some(bv) = bias_value {
        if bv <= 0.0 {
            return vec![];
        }
    }

    // Find parent of this individual
    let parent_ids = pedigree.parents_of(individual_id);
    let parent_id = match parent_ids.first() {
        Some(pid) => *pid,
        None => return vec![],
    };
    let parent_haplotype = match current_haplotypes.get(parent_id) {
        Some(h) => h,
        None => return vec![],
    };

    // Find known descendants of this individual
    let desc_ids = descendants(pedigree, individual_id);
    let known_desc: Vec<&str> = desc_ids
        .iter()
        .filter(|id| id.as_str() != individual_id)
        .filter_map(|id| pedigree.get_individual_by_id(id))
        .filter(|ind| {
            matches!(
                ind.haplotype_class,
                HaplotypeClass::Known | HaplotypeClass::Suspect | HaplotypeClass::Fixed
            )
        })
        .map(|ind| ind.id.as_str())
        .collect();

    if known_desc.is_empty() {
        return vec![];
    }

    // Compute distance to MRCA of known descendants
    let mrca_id = if known_desc.len() == 1 {
        known_desc[0].to_owned()
    } else {
        let mut current_mrca = known_desc[0].to_owned();
        for &other in &known_desc[1..] {
            if let Some(m) = most_recent_common_ancestor(pedigree, &current_mrca, other) {
                current_mrca = m;
            }
        }
        current_mrca
    };

    let mut distance_to_mrca = shortest_path_length(pedigree, individual_id, &mrca_id)
        .unwrap_or(0);
    // If the MRCA is itself unknown, add one to the distance
    if let Some(mrca_ind) = pedigree.get_individual_by_id(&mrca_id) {
        if mrca_ind.haplotype_class == HaplotypeClass::Unknown {
            distance_to_mrca += 1;
        }
    }
    let distance_to_mrca = distance_to_mrca.max(1) as u32;

    // Use SimpleDifferenceMatrix to determine mutation direction per marker
    use crate::hungarian::SimpleDifferenceMatrix;

    let mut biases = Vec::new();
    for marker in &marker_set.markers {
        let parent_alleles: Vec<_> = parent_haplotype
            .get_alleles_by_marker_name(&marker.name);

        // Collect mutation direction tuples for each known descendant
        let mutation_lists: Vec<Vec<&'static str>> = known_desc
            .iter()
            .filter_map(|&desc_id| {
                let desc_hap = current_haplotypes.get(desc_id)?;
                let desc_alleles: Vec<_> = desc_hap.get_alleles_by_marker_name(&marker.name);
                if parent_alleles.is_empty() || desc_alleles.is_empty() {
                    return None;
                }
                if parent_alleles.len() != desc_alleles.len() {
                    return None;
                }
                let matrix = SimpleDifferenceMatrix::new(parent_alleles.clone(), desc_alleles);
                let deltas = matrix.calculate_mutations();
                let dirs: Vec<&'static str> = deltas
                    .iter()
                    .map(|d| {
                        if d.step > 0 { "up" }
                        else if d.step < 0 { "down" }
                        else { "none" }
                    })
                    .collect();
                Some(dirs)
            })
            .collect();

        if mutation_lists.is_empty() {
            continue;
        }

        // All descendants must agree on the direction for each copy
        let first = &mutation_lists[0];
        let all_agree = mutation_lists.iter().all(|ml| ml == first);
        if !all_agree {
            continue;
        }

        for (copy_nr, &direction) in first.iter().enumerate() {
            if direction == "none" {
                continue;
            }
            let target_mass = match bias_value {
                Some(bv) => bv,
                None => crate::simulation::bias::default_bias_target_mass(distance_to_mrca),
            };
            biases.push(Bias {
                marker: marker.clone(),
                copy_nr: copy_nr as u32,
                direction: if direction == "up" { BiasDirection::Up } else { BiasDirection::Down },
                target_mass,
            });
        }
    }

    biases
}

// ---------------------------------------------------------------------------
// simulate_pedigree_probability_batch
// ---------------------------------------------------------------------------

/// Simulate one batch of `batch_length` iterations estimating the average
/// pedigree probability (Step 1 of the simulation).
///
/// For each iteration:
/// 1. Traverse BFS-ordered unknown individuals from root.
/// 2. For each unknown, sample a mutated haplotype from its parent using biased
///    mutation.
/// 3. Collect "fixed" edge probabilities for relationships not traversed.
/// 4. Accumulate importance-weighted sample.
pub fn simulate_pedigree_probability_batch(
    pedigree: &Pedigree,
    root_id: &str,
    marker_set: &MarkerSet,
    params: &SimulationParameters,
    rng: &mut impl Rng,
    batch_length: u64,
) -> Result<BatchResult> {
    let mut result = BatchResult::new();

    // BFS traversal order from root (flatten layers)
    let layers = bfs_layers(pedigree, root_id)?;
    let bfs_order: Vec<String> = layers.into_iter().flatten().collect();

    // child → parent map (single parent per individual in a Y-STR pedigree)
    let child_of: HashMap<String, String> = pedigree
        .relationships
        .iter()
        .map(|r| (r.child_id.clone(), r.parent_id.clone()))
        .collect();

    // Ordered unknown IDs in BFS order
    let ordered_unknown: Vec<String> = bfs_order
        .iter()
        .filter(|id| {
            pedigree
                .get_individual_by_id(id)
                .map(|ind| ind.haplotype_class == HaplotypeClass::Unknown)
                .unwrap_or(false)
        })
        .cloned()
        .collect();

    // Base haplotype map from known individuals
    let base_haplotypes: HashMap<String, Haplotype> = pedigree
        .individuals
        .iter()
        .map(|ind| (ind.id.clone(), ind.haplotype.clone()))
        .collect();

    for _ in 0..batch_length {
        let mut haplotypes = base_haplotypes.clone();
        let mut w_total = 1.0f64;
        let mut u_total = 1.0f64;
        let mut simulated_edges: std::collections::HashSet<(String, String)> =
            std::collections::HashSet::new();

        for uid in &ordered_unknown {
            let pid = match child_of.get(uid) {
                Some(p) => p,
                None => continue,
            };
            let parent_hap = match haplotypes.get(pid) {
                Some(h) => h.clone(),
                None => continue,
            };

            let biases = get_biases_for_individual(
                pedigree,
                uid,
                marker_set,
                &haplotypes,
                params.bias,
            );

            let (new_hap, w, u) = mutate_haplotype(
                &parent_hap,
                marker_set,
                params.two_step_mutation_fraction,
                &biases,
                rng,
            );

            haplotypes.insert(uid.clone(), new_hap);
            w_total *= w;
            u_total *= u;
            simulated_edges.insert((pid.clone(), uid.clone()));
        }

        // Compute edge probabilities for relationships NOT traversed by sampling
        let unused_probs: Vec<f64> = pedigree
            .relationships
            .iter()
            .filter(|r| {
                !simulated_edges.contains(&(r.parent_id.clone(), r.child_id.clone()))
                    && !simulated_edges.contains(&(r.child_id.clone(), r.parent_id.clone()))
            })
            .filter_map(|r| {
                let parent_hap = haplotypes.get(&r.parent_id)?;
                let child_hap = haplotypes.get(&r.child_id)?;
                Some(get_edge_probability(
                    parent_hap,
                    child_hap,
                    marker_set,
                    params.two_step_mutation_fraction,
                ))
            })
            .collect();

        // If any fixed edge has probability 0, skip this iteration
        if unused_probs.iter().any(|&p| p == 0.0) {
            // weight_sum still advances to avoid bias in the denominator
            let importance_weight = if w_total == 0.0 {
                Decimal::ZERO
            } else {
                Decimal::try_from(u_total / w_total).unwrap_or(Decimal::ZERO)
            };
            result.accumulate(Decimal::ZERO, importance_weight);
            continue;
        }

        let ped_prob: f64 = unused_probs.iter().product();
        let importance_weight = if w_total == 0.0 {
            Decimal::ZERO
        } else {
            Decimal::try_from(u_total / w_total).unwrap_or(Decimal::ZERO)
        };
        let probability = Decimal::try_from(ped_prob).unwrap_or(Decimal::ZERO);

        result.accumulate(probability, importance_weight);
    }

    Ok(result)
}

// ---------------------------------------------------------------------------
// simulate_matching_haplotypes_batch
// ---------------------------------------------------------------------------

/// Simulate one batch estimating match probabilities (Step 2 of the simulation).
///
/// For each iteration:
/// 1. Pick one "fixed" unknown individual weighted by picking_probabilities.
/// 2. Assign the suspect haplotype to that individual (treating it as "known").
/// 3. Simulate all other unknowns from their parents.
/// 4. Count how many unknowns match the suspect haplotype.
/// 5. Accumulate the conditional probability.
pub fn simulate_matching_haplotypes_batch(
    pedigree: &Pedigree,
    root_id: &str,
    suspect_haplotype: &Haplotype,
    avg_pedigree_probability: Decimal,
    marker_set: &MarkerSet,
    params: &SimulationParameters,
    is_outside: bool,
    rng: &mut impl Rng,
    batch_length: u64,
) -> Result<BatchResult> {
    let mut result = BatchResult::new();

    // Build picking probability distribution
    let unknown_ids: Vec<String> = pedigree
        .individuals
        .iter()
        .filter(|i| i.haplotype_class == HaplotypeClass::Unknown)
        .map(|i| i.id.clone())
        .collect();

    let picking_probs: Vec<f64> = unknown_ids
        .iter()
        .map(|id| {
            pedigree
                .picking_probabilities
                .get(id)
                .and_then(|d| f64::try_from(*d).ok())
                .unwrap_or(0.0)
        })
        .collect();

    if picking_probs.iter().all(|&p| p == 0.0) || unknown_ids.is_empty() {
        return Ok(result);
    }

    use rand::distributions::WeightedIndex;
    let pick_dist = match WeightedIndex::new(&picking_probs) {
        Ok(d) => d,
        Err(_) => return Ok(result),
    };

    // BFS traversal from root
    let layers = bfs_layers(pedigree, root_id)?;
    let bfs_order: Vec<String> = layers.into_iter().flatten().collect();

    let ordered_unknown: Vec<String> = bfs_order
        .iter()
        .filter(|id| unknown_ids.contains(id))
        .cloned()
        .collect();

    let child_of: HashMap<String, String> = pedigree
        .relationships
        .iter()
        .map(|r| (r.child_id.clone(), r.parent_id.clone()))
        .collect();

    let base_haplotypes: HashMap<String, Haplotype> = pedigree
        .individuals
        .iter()
        .map(|ind| (ind.id.clone(), ind.haplotype.clone()))
        .collect();

    for _ in 0..batch_length {
        // Pick the fixed individual
        let fixed_idx = pick_dist.sample(rng);
        let fixed_id = &unknown_ids[fixed_idx];

        // Normalise picking probability
        let total_pick: f64 = picking_probs.iter().sum();
        let fixed_prob = picking_probs[fixed_idx] / total_pick;

        // Build haplotype map: fixed individual gets suspect haplotype
        let mut haplotypes = base_haplotypes.clone();
        haplotypes.insert(fixed_id.clone(), suspect_haplotype.clone());

        let mut w_total = 1.0f64;
        let mut u_total = 1.0f64;
        let mut simulated_edges: std::collections::HashSet<(String, String)> =
            std::collections::HashSet::new();

        // Simulate all unknown individuals except the fixed one
        for uid in &ordered_unknown {
            if uid == fixed_id {
                continue;
            }

            let pid = match child_of.get(uid) {
                Some(p) => p,
                None => continue,
            };
            let parent_hap = match haplotypes.get(pid) {
                Some(h) => h.clone(),
                None => continue,
            };

            let biases = get_biases_for_individual(
                pedigree,
                uid,
                marker_set,
                &haplotypes,
                params.bias,
            );

            let (new_hap, w, u) = mutate_haplotype(
                &parent_hap,
                marker_set,
                params.two_step_mutation_fraction,
                &biases,
                rng,
            );

            haplotypes.insert(uid.clone(), new_hap);
            w_total *= w;
            u_total *= u;
            simulated_edges.insert((pid.clone(), uid.clone()));
        }

        // All edge probabilities (all relationships)
        let all_probs: Vec<f64> = pedigree
            .relationships
            .iter()
            .filter_map(|r| {
                let parent_hap = haplotypes.get(&r.parent_id)?;
                let child_hap = haplotypes.get(&r.child_id)?;
                Some(get_edge_probability(
                    parent_hap,
                    child_hap,
                    marker_set,
                    params.two_step_mutation_fraction,
                ))
            })
            .collect();

        let ped_prob: f64 = all_probs.iter().product();

        // Which simulated edges contribute to the simulation probability?
        let sim_prob_factor: f64 = pedigree
            .relationships
            .iter()
            .filter(|r| simulated_edges.contains(&(r.parent_id.clone(), r.child_id.clone())))
            .filter_map(|r| {
                let parent_hap = haplotypes.get(&r.parent_id)?;
                let child_hap = haplotypes.get(&r.child_id)?;
                Some(get_edge_probability(
                    parent_hap,
                    child_hap,
                    marker_set,
                    params.two_step_mutation_fraction,
                ))
            })
            .product();

        let simulation_probability = fixed_prob * sim_prob_factor;

        // Conditional probability
        let cond: f64 = if avg_pedigree_probability == Decimal::ZERO || avg_pedigree_probability == Decimal::ZERO {
            0.0
        } else {
            let avg = f64::try_from(avg_pedigree_probability).unwrap_or(0.0);
            if avg == 0.0 { 0.0 } else { ped_prob / avg }
        };

        // Count matches (non-excluded unknowns that equal suspect haplotype)
        let non_excl_matches: Vec<String> = ordered_unknown
            .iter()
            .filter(|uid| {
                let matches_hap = haplotypes
                    .get(uid.as_str())
                    .map(|h| h == suspect_haplotype)
                    .unwrap_or(false);
                let not_excluded = pedigree
                    .get_individual_by_id(uid)
                    .map(|ind| !ind.exclude)
                    .unwrap_or(false);
                matches_hap && not_excluded
            })
            .cloned()
            .collect();

        let total_matching: std::collections::HashSet<String> = {
            let mut s: std::collections::HashSet<String> = ordered_unknown
                .iter()
                .filter(|uid| {
                    haplotypes
                        .get(uid.as_str())
                        .map(|h| h == suspect_haplotype)
                        .unwrap_or(false)
                })
                .cloned()
                .collect();
            s.insert(fixed_id.clone());
            s
        };
        let total_number = total_matching.len();

        let probability = if is_outside {
            if simulation_probability == 0.0 {
                0.0
            } else {
                cond / simulation_probability
            }
        } else {
            let n_matches = non_excl_matches.len();
            if n_matches >= 1 && simulation_probability != 0.0 && total_number > 0 {
                cond / (simulation_probability * total_number as f64)
            } else {
                0.0
            }
        };

        let importance_weight = if w_total == 0.0 {
            Decimal::ZERO
        } else {
            Decimal::try_from(u_total / w_total).unwrap_or(Decimal::ZERO)
        };

        let probability_dec = Decimal::try_from(probability).unwrap_or(Decimal::ZERO);

        // Accumulate match-count data
        if !is_outside && !non_excl_matches.is_empty() {
            let k = non_excl_matches.len() as u32;
            *result.match_accumulators.entry(k).or_default() +=
                probability_dec * importance_weight;
        }

        // Accumulate per-individual marginals
        for uid in &non_excl_matches {
            *result.per_individual.entry(uid.clone()).or_default() +=
                probability_dec * importance_weight;
        }

        result.accumulate(probability_dec, importance_weight);
    }

    Ok(result)
}
