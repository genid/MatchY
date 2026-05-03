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
use crate::simulation::mutation::{mutate_haplotype_precomputed, neutral_step_probabilities};
use crate::simulation::probability::get_edge_probability;
use crate::graph::{bfs_layers, descendants, most_recent_common_ancestor, shortest_path_length};
use crate::Result;
use rand::prelude::*;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use rust_decimal::Decimal;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// BatchResult
// ---------------------------------------------------------------------------

/// Running statistics for one Monte Carlo batch (many iterations).
///
/// `weighted_sum` and `weight_sum` are kept as `f64` so that large IS weights
/// (which can exceed `rust_decimal::Decimal`'s ~7.9×10²⁸ ceiling when a fixed
/// bias is applied uniformly across many unknowns) are accumulated without
/// silent overflow.  The running *mean* (a probability in `[0,1]`) is only
/// converted to `Decimal` at the point it is appended to `running_means`.
#[derive(Debug, Clone)]
pub struct BatchResult {
    /// Sum of (probability × importance_weight) — numerator, kept as f64
    pub weighted_sum: f64,
    /// Sum of importance_weights — denominator, kept as f64
    pub weight_sum: f64,
    /// Number of iterations completed
    pub iterations: u64,
    /// Per-iteration running means — mirrors Python's model_probabilities[m] list.
    /// Used by the convergence check, which requires ALL running means to be within
    /// threshold (not just the final one).  Seeded with the last mean of the previous
    /// trial so that the Python `model_probabilities[m][-1:]` carryover is reproduced.
    pub running_means: Vec<Decimal>,
    /// Per-match-count weighted accumulators: match_count → weighted_sum (f64)
    pub match_accumulators: HashMap<u32, f64>,
    /// Per-individual weighted accumulators: individual_id → weighted_sum (f64)
    pub per_individual: HashMap<String, f64>,
}

impl BatchResult {
    pub fn new() -> Self {
        Self {
            weighted_sum: 0.0,
            weight_sum: 0.0,
            iterations: 0,
            running_means: Vec::new(),
            match_accumulators: HashMap::new(),
            per_individual: HashMap::new(),
        }
    }

    /// Create a BatchResult pre-seeded with a carryover mean from the previous trial.
    /// Mirrors Python: `model_probabilities = {m: model_probabilities[m][-1:] for m in range(3)}`.
    pub fn new_with_seed(seed_mean: Decimal) -> Self {
        Self {
            running_means: vec![seed_mean],
            ..Self::new()
        }
    }

    /// Running importance-weighted mean estimate.
    /// Returns `None` when no iterations have been accumulated yet.
    pub fn running_mean(&self) -> Option<Decimal> {
        if self.weight_sum == 0.0 {
            return None;
        }
        let mean = self.weighted_sum / self.weight_sum;
        Decimal::try_from(mean).ok()
    }

    /// Add one iteration's result, appending the running mean to the history.
    pub fn accumulate(&mut self, probability: f64, importance_weight: f64) {
        self.weight_sum += importance_weight;
        self.weighted_sum += probability * importance_weight;
        self.iterations += 1;
        if let Some(mean) = self.running_mean() {
            self.running_means.push(mean);
        }
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
// fixed_id: treat this individual as "Known" during bias computation — mirrors
// Python simulation.py:439-441 where fixed_individual_id.haplotype_class is
// temporarily set to "known" so that bias logic sees it as a known descendant.
// Overlay lookup: simulated unknowns shadow base haplotypes.
#[inline]
fn hap_lookup<'a>(
    base: &'a HashMap<String, Haplotype>,
    sim: &'a HashMap<String, Haplotype>,
    id: &str,
) -> Option<&'a Haplotype> {
    sim.get(id).or_else(|| base.get(id))
}

fn get_biases_for_individual(
    pedigree: &Pedigree,
    individual_id: &str,
    marker_set: &MarkerSet,
    base_haps: &HashMap<String, Haplotype>,
    sim_haps: &HashMap<String, Haplotype>,
    bias_value: Option<f64>,
    fixed_id: Option<&str>,
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
    let parent_haplotype = match hap_lookup(base_haps, sim_haps, parent_id) {
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
            let class_known = matches!(
                ind.haplotype_class,
                HaplotypeClass::Known | HaplotypeClass::Suspect | HaplotypeClass::Fixed
            );
            // Treat fixed individual as "known" — mirrors Python simulation.py:439-441
            let is_fixed_override = fixed_id.map(|fid| fid == ind.id.as_str()).unwrap_or(false);
            class_known || is_fixed_override
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
                let desc_hap = hap_lookup(base_haps, sim_haps, desc_id)?;
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
    initial_sums: Option<(f64, f64)>,
) -> Result<BatchResult> {
    let mut result = if let Some((ws, w)) = initial_sums {
        BatchResult { weighted_sum: ws, weight_sum: w, ..BatchResult::new() }
    } else {
        BatchResult::new()
    };

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

    // (C) Precompute neutral step probability tables — constant for this batch
    let neutral_tables: HashMap<String, [f64; 5]> = marker_set
        .markers
        .iter()
        .map(|m| (m.name.clone(), neutral_step_probabilities(m.single_copy_mutation_rate(), params.two_step_mutation_fraction)))
        .collect();

    // (D) Unknown ID set for classifying pedigree edges
    let unknown_id_set: std::collections::HashSet<&str> =
        ordered_unknown.iter().map(|s| s.as_str()).collect();

    // (D+5) Precompute known→known edge log-probability sum — constant across all iterations.
    // Log-space prevents underflow for large pedigrees with many required mutations.
    let log_const_edge_prob: f64 = pedigree
        .relationships
        .iter()
        .filter(|r| {
            !unknown_id_set.contains(r.parent_id.as_str())
                && !unknown_id_set.contains(r.child_id.as_str())
        })
        .filter_map(|r| {
            let ph = base_haplotypes.get(&r.parent_id)?;
            let ch = base_haplotypes.get(&r.child_id)?;
            Some(get_edge_probability(ph, ch, marker_set, params.two_step_mutation_fraction).ln())
        })
        .sum();

    // (D) Prefilter relationships where parent is simulated-unknown and child is known —
    // these must be evaluated per-iteration since the parent haplotype varies.
    let var_edge_rels: Vec<(&str, &str)> = pedigree
        .relationships
        .iter()
        .filter(|r| {
            unknown_id_set.contains(r.parent_id.as_str())
                && !unknown_id_set.contains(r.child_id.as_str())
        })
        .map(|r| (r.parent_id.as_str(), r.child_id.as_str()))
        .collect();

    // One RNG per worker thread — avoids creating batch_length separate RNG objects
    let n_chunks = rayon::current_num_threads().max(1);
    let chunk_size = (batch_length as usize + n_chunks - 1) / n_chunks;
    let chunk_seeds: Vec<u64> = (0..n_chunks).map(|_| rng.gen()).collect();

    let all_samples: Vec<Vec<(f64, f64)>> = chunk_seeds
        .into_par_iter()
        .enumerate()
        .map(|(chunk_idx, seed)| {
            let mut local_rng = SmallRng::seed_from_u64(seed);
            let start = chunk_idx * chunk_size;
            let end = ((chunk_idx + 1) * chunk_size).min(batch_length as usize);
            (start..end)
                .map(|_| {
                    // (B) Overlay: only allocate storage for simulated unknowns; look up
                    // known individuals directly from the immutable base map.
                    let mut sim: HashMap<String, Haplotype> =
                        HashMap::with_capacity(ordered_unknown.len());
                    let mut w_total = 1.0f64;
                    let mut u_total = 1.0f64;

                    for uid in &ordered_unknown {
                        let pid = match child_of.get(uid) {
                            Some(p) => p,
                            None => continue,
                        };
                        let parent_hap = match hap_lookup(&base_haplotypes, &sim, pid) {
                            Some(h) => h.clone(),
                            None => continue,
                        };

                        let biases = get_biases_for_individual(
                            pedigree, uid, marker_set, &base_haplotypes, &sim, params.bias, None,
                        );

                        // (C) Use precomputed neutral step tables
                        let (new_hap, w, u) = mutate_haplotype_precomputed(
                            &parent_hap, marker_set, &neutral_tables, &biases, &mut local_rng,
                        );

                        sim.insert(uid.clone(), new_hap);
                        w_total *= w;
                        u_total *= u;
                    }

                    let importance_weight = if w_total == 0.0 {
                        0.0f64
                    } else {
                        let iw = u_total / w_total;
                        if iw.is_finite() { iw } else { 0.0 }
                    };

                    // (D+5) Known→known edges: -inf means at least one probability is zero.
                    if !log_const_edge_prob.is_finite() {
                        return (0.0f64, importance_weight);
                    }

                    // (D+5) Variable edges in log-space — prevents underflow for deep pedigrees.
                    let log_var_edge: f64 = var_edge_rels
                        .iter()
                        .filter_map(|&(pid, cid)| {
                            let ph = sim.get(pid)?;
                            let ch = base_haplotypes.get(cid)?;
                            Some(get_edge_probability(ph, ch, marker_set, params.two_step_mutation_fraction).ln())
                        })
                        .sum();

                    if !log_var_edge.is_finite() {
                        return (0.0f64, importance_weight);
                    }

                    let ped_prob = (log_const_edge_prob + log_var_edge).exp();
                    let probability = if ped_prob.is_finite() { ped_prob } else { 0.0 };
                    (probability, importance_weight)
                })
                .collect()
        })
        .collect();

    for chunk in all_samples {
        for (probability, importance_weight) in chunk {
            result.accumulate(probability, importance_weight);
        }
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
    // Full carry from previous trial — all accumulators are cumulative across trials,
    // matching Python where weight_sums/weighted_sums/per_individual_weighted_sums
    // are never reset between trials.
    initial_carry: Option<BatchResult>,
) -> Result<BatchResult> {
    let mut result = if let Some(mut carry) = initial_carry {
        // Reset per-trial fields; cumulative accumulators are preserved.
        carry.running_means.clear();
        carry.iterations = 0;
        carry
    } else {
        BatchResult::new()
    };

    // All unknowns are candidates to be "fixed" (matching Python's picking pool behaviour).
    // The exclude flag only affects match counting, not the picking pool.
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

    // (C) Precompute neutral step probability tables — constant for this batch
    let neutral_tables_match: HashMap<String, [f64; 5]> = marker_set
        .markers
        .iter()
        .map(|m| (m.name.clone(), neutral_step_probabilities(m.single_copy_mutation_rate(), params.two_step_mutation_fraction)))
        .collect();

    // (D) Unknown ID set for classifying pedigree edges
    let unknown_id_set_match: std::collections::HashSet<&str> =
        unknown_ids.iter().map(|s| s.as_str()).collect();

    // (D+5) Precompute known→known edge log-probability sum — constant across all iterations.
    // Log-space prevents underflow for large pedigrees with many required mutations.
    let log_const_edge_prob_match: f64 = pedigree
        .relationships
        .iter()
        .filter(|r| {
            !unknown_id_set_match.contains(r.parent_id.as_str())
                && !unknown_id_set_match.contains(r.child_id.as_str())
        })
        .filter_map(|r| {
            let ph = base_haplotypes.get(&r.parent_id)?;
            let ch = base_haplotypes.get(&r.child_id)?;
            Some(get_edge_probability(ph, ch, marker_set, params.two_step_mutation_fraction).ln())
        })
        .sum();

    let total_pick: f64 = picking_probs.iter().sum();
    let avg_pp = f64::try_from(avg_pedigree_probability).unwrap_or(0.0);

    // One RNG per worker thread — avoids creating batch_length separate RNG objects
    let n_chunks = rayon::current_num_threads().max(1);
    let chunk_size = (batch_length as usize + n_chunks - 1) / n_chunks;
    let chunk_seeds: Vec<u64> = (0..n_chunks).map(|_| rng.gen()).collect();

    // Each item: (probability, importance_weight, match_acc_deltas, per_individual_deltas)
    type IterSample = (f64, f64, Vec<(u32, f64)>, Vec<(String, f64)>);

    let all_samples: Vec<Vec<IterSample>> = chunk_seeds
        .into_par_iter()
        .enumerate()
        .map(|(chunk_idx, seed)| {
            let mut local_rng = SmallRng::seed_from_u64(seed);
            let start = chunk_idx * chunk_size;
            let end = ((chunk_idx + 1) * chunk_size).min(batch_length as usize);
            (start..end)
                .map(|_| {
                    let fixed_idx = pick_dist.sample(&mut local_rng);
                    let fixed_id = &unknown_ids[fixed_idx];
                    let fixed_prob = picking_probs[fixed_idx] / total_pick;

                    // (B) Overlay: assign suspect haplotype to fixed_id; all other unknowns
                    // are added incrementally during BFS. Known individuals are read directly
                    // from the immutable base map — no full clone needed.
                    let mut sim: HashMap<String, Haplotype> =
                        HashMap::with_capacity(ordered_unknown.len());
                    sim.insert(fixed_id.clone(), suspect_haplotype.clone());

                    let mut w_total = 1.0f64;
                    let mut u_total = 1.0f64;
                    let mut simulated_edges: std::collections::HashSet<(String, String)> =
                        std::collections::HashSet::new();

                    for uid in &ordered_unknown {
                        if uid == fixed_id {
                            continue;
                        }
                        let pid = match child_of.get(uid) {
                            Some(p) => p,
                            None => continue,
                        };
                        let parent_hap = match hap_lookup(&base_haplotypes, &sim, pid) {
                            Some(h) => h.clone(),
                            None => continue,
                        };
                        let biases = get_biases_for_individual(
                            pedigree, uid, marker_set, &base_haplotypes, &sim, params.bias, Some(fixed_id.as_str()),
                        );
                        // (C) Use precomputed neutral step tables
                        let (new_hap, w, u) = mutate_haplotype_precomputed(
                            &parent_hap, marker_set, &neutral_tables_match, &biases, &mut local_rng,
                        );
                        sim.insert(uid.clone(), new_hap);
                        w_total *= w;
                        u_total *= u;
                        simulated_edges.insert((pid.clone(), uid.clone()));
                    }

                    // (E + D + 5) Single-pass edge probability computation in log-space.
                    // Start log_pp from the precomputed known→known log-sum; skip those
                    // edges in the loop to avoid recomputing them.
                    let (ped_prob, sim_prob_factor) = {
                        let mut log_pp = log_const_edge_prob_match;
                        let mut log_sp = 0.0f64;
                        for r in &pedigree.relationships {
                            // (D) Skip known→known — already in log_const_edge_prob_match
                            if !unknown_id_set_match.contains(r.parent_id.as_str())
                                && !unknown_id_set_match.contains(r.child_id.as_str())
                            {
                                continue;
                            }
                            let ph = match hap_lookup(&base_haplotypes, &sim, &r.parent_id) { Some(h) => h, None => continue };
                            let ch = match hap_lookup(&base_haplotypes, &sim, &r.child_id) { Some(h) => h, None => continue };
                            let ln_p = get_edge_probability(ph, ch, marker_set, params.two_step_mutation_fraction).ln();
                            log_pp += ln_p;
                            // (E) Accumulate simulated-edge factor in the same pass
                            if simulated_edges.contains(&(r.parent_id.clone(), r.child_id.clone())) {
                                log_sp += ln_p;
                            }
                        }
                        (log_pp.exp(), log_sp.exp())
                    };

                    let simulation_probability = fixed_prob * sim_prob_factor;
                    let cond = if avg_pp == 0.0 { 0.0 } else { ped_prob / avg_pp };

                    let non_excl_matches: Vec<String> = ordered_unknown
                        .iter()
                        .filter(|uid| {
                            hap_lookup(&base_haplotypes, &sim, uid.as_str()).map(|h| h == suspect_haplotype).unwrap_or(false)
                                && pedigree.get_individual_by_id(uid).map(|i| !i.exclude).unwrap_or(false)
                        })
                        .cloned()
                        .collect();

                    let total_number = {
                        let mut s: std::collections::HashSet<&str> = ordered_unknown
                            .iter()
                            .filter(|uid| hap_lookup(&base_haplotypes, &sim, uid.as_str()).map(|h| h == suspect_haplotype).unwrap_or(false))
                            .map(|s| s.as_str())
                            .collect();
                        s.insert(fixed_id.as_str());
                        s.len()
                    };

                    let probability = if is_outside {
                        if simulation_probability == 0.0 { 0.0 } else { cond / simulation_probability }
                    } else {
                        let n = non_excl_matches.len();
                        if n >= 1 && simulation_probability != 0.0 && total_number > 0 {
                            cond / (simulation_probability * total_number as f64)
                        } else {
                            0.0
                        }
                    };

                    let importance_weight = if w_total == 0.0 {
                        0.0f64
                    } else {
                        let iw = u_total / w_total;
                        if iw.is_finite() { iw } else { 0.0 }
                    };
                    let probability = if probability.is_finite() { probability } else { 0.0 };
                    let weighted = probability * importance_weight;

                    let mut match_acc: Vec<(u32, f64)> = Vec::new();
                    if !is_outside && !non_excl_matches.is_empty() {
                        match_acc.push((non_excl_matches.len() as u32, weighted));
                    }
                    let per_ind: Vec<(String, f64)> = non_excl_matches
                        .into_iter()
                        .map(|uid| (uid, weighted))
                        .collect();

                    (probability, importance_weight, match_acc, per_ind)
                })
                .collect()
        })
        .collect();

    for chunk in all_samples {
        for (prob, weight, match_acc, per_ind) in chunk {
            for (k, v) in match_acc {
                *result.match_accumulators.entry(k).or_default() += v;
            }
            for (id, v) in per_ind {
                *result.per_individual.entry(id).or_default() += v;
            }
            result.accumulate(prob, weight);
        }
    }

    Ok(result)
}
