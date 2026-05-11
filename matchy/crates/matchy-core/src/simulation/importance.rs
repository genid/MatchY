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
use crate::simulation::probability::{calculate_mutation_probability, get_edge_probability};
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
/// `weighted_sum` and `weight_sum` are kept as `f64` to avoid overflow or precision
/// loss. `running_means` also stores `f64` — using `Decimal` would silently drop
/// probabilities smaller than ~1e-28 (Decimal's representable floor), causing the
/// convergence check to stall forever on large pedigrees.
#[derive(Debug, Clone)]
pub struct BatchResult {
    /// Sum of (probability × importance_weight) — numerator, kept as f64
    pub weighted_sum: f64,
    /// Sum of importance_weights — denominator, kept as f64
    pub weight_sum: f64,
    /// Sum of importance_weight² — for effective sample size computation
    pub sum_iw_sq: f64,
    /// Number of iterations completed
    pub iterations: u64,
    /// Per-iteration running means stored as f64 — mirrors Python's model_probabilities[m].
    /// Used by the convergence check. Seeded with the last mean of the previous trial.
    pub running_means: Vec<f64>,
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
            sum_iw_sq: 0.0,
            iterations: 0,
            running_means: Vec::new(),
            match_accumulators: HashMap::new(),
            per_individual: HashMap::new(),
        }
    }

    /// Create a BatchResult pre-seeded with a carryover mean from the previous trial.
    /// Mirrors Python: `model_probabilities = {m: model_probabilities[m][-1:] for m in range(3)}`.
    pub fn new_with_seed(seed_mean: f64) -> Self {
        Self {
            running_means: vec![seed_mean],
            ..Self::new()
        }
    }

    /// Running importance-weighted mean estimate as f64.
    /// Returns 0.0 when no iterations have been accumulated yet.
    pub fn running_mean_f64(&self) -> f64 {
        if self.weight_sum == 0.0 { 0.0 } else { self.weighted_sum / self.weight_sum }
    }

    /// Running mean as Decimal for display/reporting. Returns None for zero or
    /// sub-Decimal-precision values (below ~1e-28).
    pub fn running_mean(&self) -> Option<Decimal> {
        if self.weight_sum == 0.0 { return None; }
        Decimal::try_from(self.weighted_sum / self.weight_sum).ok()
    }

    /// Add one iteration's result, appending the running mean to the history.
    pub fn accumulate(&mut self, probability: f64, importance_weight: f64) {
        self.weight_sum += importance_weight;
        self.weighted_sum += probability * importance_weight;
        self.sum_iw_sq += importance_weight * importance_weight;
        self.iterations += 1;
        if self.weight_sum > 0.0 {
            self.running_means.push(self.weighted_sum / self.weight_sum);
        }
    }

    /// Effective sample size: weight_sum² / sum_iw_sq.
    /// A low ESS (relative to iterations) indicates IS degeneracy.
    pub fn effective_sample_size(&self) -> f64 {
        if self.sum_iw_sq == 0.0 { return 0.0; }
        self.weight_sum * self.weight_sum / self.sum_iw_sq
    }

    /// ESS as a fraction of total iterations in [0, 1]. Values < 0.01 indicate severe degeneracy.
    pub fn relative_ess(&self) -> f64 {
        if self.iterations == 0 { return 0.0; }
        self.effective_sample_size() / self.iterations as f64
    }

    /// True when all non-zero importance weights correspond to zero-probability samples —
    /// the clearest sign of IS degeneracy (numerator is 0 but denominator is not).
    pub fn is_degenerate(&self) -> bool {
        self.weight_sum > 0.0 && self.weighted_sum == 0.0 && self.iterations > 0
    }
}

impl Default for BatchResult {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Zero-probability debug capture
// ---------------------------------------------------------------------------

struct ZeroProbMarkerInfo {
    name: String,
    parent_vals: Vec<String>,
    child_vals: Vec<String>,
    probability: f64,
}

struct ZeroProbEdgeInfo {
    parent_id: String,
    child_id: String,
    markers: Vec<ZeroProbMarkerInfo>,
}

struct ZeroProbUnknownHap {
    id: String,
    parent_id: String,
    /// (marker_name, [(allele_display, step_taken)])
    markers: Vec<(String, Vec<(String, i32)>)>,
}

struct ZeroProbSample {
    unknowns: Vec<ZeroProbUnknownHap>,
    /// Only the edges whose overall probability is 0
    zero_edges: Vec<ZeroProbEdgeInfo>,
}

/// Snapshot the current simulation state for one zero-probability iteration.
fn collect_zero_prob_sample(
    sim: &HashMap<String, Haplotype>,
    base_haplotypes: &HashMap<String, Haplotype>,
    ordered_unknown: &[String],
    child_of: &HashMap<String, String>,
    var_edge_rels: &[(&str, &str)],
    marker_set: &MarkerSet,
    two_step_fraction: f64,
) -> ZeroProbSample {
    let unknowns = ordered_unknown
        .iter()
        .filter_map(|uid| {
            let hap = sim.get(uid.as_str())?;
            let parent_id = child_of.get(uid.as_str()).cloned().unwrap_or_default();
            let markers = marker_set
                .markers
                .iter()
                .filter_map(|m| {
                    let alleles = hap.alleles.get(&m.name)?;
                    let vals: Vec<(String, i32)> = alleles
                        .iter()
                        .map(|a| (a.display(), a.mutation_value.unwrap_or(0)))
                        .collect();
                    Some((m.name.clone(), vals))
                })
                .collect();
            Some(ZeroProbUnknownHap { id: uid.clone(), parent_id, markers })
        })
        .collect();

    let zero_edges = var_edge_rels
        .iter()
        .filter_map(|&(pid, cid)| {
            let ph = sim.get(pid)?;
            let ch = base_haplotypes.get(cid)?;
            // Only include this edge if its overall probability is 0
            let overall = get_edge_probability(ph, ch, marker_set, two_step_fraction);
            if overall > 0.0 {
                return None;
            }
            let markers = marker_set
                .markers
                .iter()
                .filter_map(|m| {
                    let parent_alleles = ph.get_alleles_by_marker_name(&m.name);
                    let child_alleles = ch.get_alleles_by_marker_name(&m.name);
                    if parent_alleles.is_empty() || child_alleles.is_empty() {
                        return None;
                    }
                    let parent_vals = parent_alleles.iter().map(|a| a.display()).collect();
                    let child_vals = child_alleles.iter().map(|a| a.display()).collect();
                    if parent_alleles.len() != child_alleles.len() {
                        return Some(ZeroProbMarkerInfo {
                            name: m.name.clone(),
                            parent_vals,
                            child_vals,
                            probability: 0.0,
                        });
                    }
                    let probability = calculate_mutation_probability(
                        &parent_alleles,
                        &child_alleles,
                        m.single_copy_mutation_rate(),
                        two_step_fraction,
                    );
                    Some(ZeroProbMarkerInfo { name: m.name.clone(), parent_vals, child_vals, probability })
                })
                .collect();
            Some(ZeroProbEdgeInfo {
                parent_id: pid.to_string(),
                child_id: cid.to_string(),
                markers,
            })
        })
        .collect();

    ZeroProbSample { unknowns, zero_edges }
}

/// Write captured zero-probability samples to a human-readable text file.
fn write_zero_prob_debug(
    samples: &[ZeroProbSample],
    path: &std::path::Path,
) -> std::io::Result<()> {
    use std::fmt::Write as FmtWrite;
    let mut out = String::new();
    let _ = writeln!(out, "MatchY — Zero-Probability Debug Samples");
    let _ = writeln!(
        out,
        "Captured {} iteration(s) where P(haplotype configuration) = 0\n",
        samples.len()
    );
    for (i, sample) in samples.iter().enumerate() {
        let _ = writeln!(out, "=== Sample {} ===", i + 1);
        let _ = writeln!(out, "\nSimulated unknown haplotypes:");
        for unk in &sample.unknowns {
            let _ = writeln!(out, "  {} (parent: {}):", unk.id, unk.parent_id);
            for (mname, alleles) in &unk.markers {
                let parts: Vec<String> = alleles
                    .iter()
                    .map(|(v, step)| {
                        if *step == 0 {
                            v.clone()
                        } else {
                            format!("{} (step {:+})", v, step)
                        }
                    })
                    .collect();
                let _ = writeln!(out, "    {:20} {}", mname, parts.join(", "));
            }
        }
        let _ = writeln!(out, "\nZero-probability edges (unknown parent → known child):");
        for edge in &sample.zero_edges {
            let _ = writeln!(out, "  {} → {}:", edge.parent_id, edge.child_id);
            for m in &edge.markers {
                let pstr = &m.parent_vals;
                let cstr = &m.child_vals;
                let flag = if m.probability == 0.0 { "  ← ZERO" } else { "" };
                let _ = writeln!(
                    out,
                    "    {:20} parent=[{}]  child=[{}]  p={:.4e}{}",
                    m.name,
                    pstr.join(", "),
                    cstr.join(", "),
                    m.probability,
                    flag
                );
            }
        }
        let _ = writeln!(out);
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    std::fs::write(path, out)
}

// ---------------------------------------------------------------------------
// Global copy labeling
// ---------------------------------------------------------------------------

/// Pre-computed consistent copy labels for all known individuals.
///
/// Labels are assigned relative to the first known individual encountered in
/// BFS order (the "root known").  For each subsequent known individual we run
/// the Hungarian algorithm between the nearest known ancestor (already
/// labelled) and the individual, propagating labels through the optimal
/// allele pairing.
///
/// Result: `assignments[(ind_id, marker_name)][sorted_idx]` = global label
/// (0-based integer). Two known individuals that share a label at a given
/// sorted index are considered to carry the "same" ancestral copy.
pub struct CopyLabeling {
    assignments: HashMap<(String, String), Vec<usize>>,
}

impl CopyLabeling {
    pub fn compute(pedigree: &Pedigree, marker_set: &crate::MarkerSet) -> Self {
        use crate::hungarian::SimpleDifferenceMatrix;

        let mut assignments: HashMap<(String, String), Vec<usize>> = HashMap::new();

        let roots = pedigree.roots();
        let root_id = match roots.first() { Some(r) => r.to_string(), None => return Self { assignments } };
        let layers = match crate::graph::bfs_layers(pedigree, &root_id) { Ok(l) => l, Err(_) => return Self { assignments } };
        let bfs_order: Vec<String> = layers.into_iter().flatten().collect();

        for marker in &marker_set.markers {
            let mname = &marker.name;

            // Locate the root known individual with alleles for this marker
            let root_known = bfs_order.iter()
                .filter_map(|id| pedigree.get_individual_by_id(id))
                .find(|ind| {
                    matches!(ind.haplotype_class,
                        crate::HaplotypeClass::Known | crate::HaplotypeClass::Suspect | crate::HaplotypeClass::Fixed)
                    && !ind.haplotype.get_alleles_by_marker_name(mname).is_empty()
                });
            let root_known = match root_known { Some(r) => r, None => continue };

            let n = root_known.haplotype.get_alleles_by_marker_name(mname).len();
            // Identity labels for the root: 0, 1, …, n-1
            assignments.insert((root_known.id.clone(), mname.clone()), (0..n).collect());

            // Process all other known individuals in BFS order
            for id in &bfs_order {
                if id == &root_known.id { continue; }
                let ind = match pedigree.get_individual_by_id(id) { Some(i) => i, None => continue };
                if !matches!(ind.haplotype_class,
                    crate::HaplotypeClass::Known | crate::HaplotypeClass::Suspect | crate::HaplotypeClass::Fixed)
                {
                    continue;
                }
                let ind_alleles = ind.haplotype.get_alleles_by_marker_name(mname);
                if ind_alleles.len() != n { continue; }

                // Walk up the pedigree to find the nearest already-labelled ancestor
                let mut cur = id.clone();
                let mut ancestor: Option<(String, Vec<usize>)> = None;
                for _ in 0..64 {
                    let parents = pedigree.parents_of(&cur);
                    let pid = match parents.first() { Some(p) => (*p).to_string(), None => break };
                    if let Some(labels) = assignments.get(&(pid.clone(), mname.clone())) {
                        ancestor = Some((pid, labels.clone()));
                        break;
                    }
                    cur = pid;
                }
                let (anc_id, anc_labels) = match ancestor {
                    Some(a) => a,
                    None => (root_known.id.clone(), (0..n).collect()),
                };

                let anc_alleles = match pedigree.get_individual_by_id(&anc_id) {
                    Some(a) => a.haplotype.get_alleles_by_marker_name(mname),
                    None => continue,
                };
                if anc_alleles.len() != n { continue; }

                let matrix = SimpleDifferenceMatrix::new(anc_alleles, ind_alleles);
                let deltas = matrix.calculate_mutations();

                // Propagate: anc's copy at sorted_idx parent_copy (with label L)
                // → ind's copy at sorted_idx child_copy gets label L
                let mut labels = vec![usize::MAX; n];
                for delta in &deltas {
                    let anc_label = anc_labels.get(delta.parent_copy as usize)
                        .copied().unwrap_or(delta.parent_copy as usize);
                    if (delta.child_copy as usize) < n {
                        labels[delta.child_copy as usize] = anc_label;
                    }
                }
                // Fill any gaps (should not happen with a valid square assignment)
                for i in 0..n { if labels[i] == usize::MAX { labels[i] = i; } }

                assignments.insert((ind.id.clone(), mname.clone()), labels);
            }
        }

        Self { assignments }
    }

    fn get_labels(&self, individual_id: &str, marker_name: &str) -> Option<&Vec<usize>> {
        self.assignments.get(&(individual_id.to_string(), marker_name.to_string()))
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
    auto_bias_strength: f64,
    auto_bias_min: f64,
    auto_bias_max: f64,
    fixed_id: Option<&str>,
    labeling: &CopyLabeling,
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
        if parent_alleles.is_empty() { continue; }
        let n_copies = parent_alleles.len();

        // Per-label unanimity using globally pre-computed copy labels.
        //
        // For each known descendant D, we run Hungarian(parent, D) to get the
        // parent→child copy mapping. We then look up D's global label for each
        // of its copies and group by label.  A label gets a bias only when ALL
        // known descendants agree on that label's mutation direction.
        //
        // This is more robust than raw per-copy unanimity: when copy values have
        // drifted in sorted order across branches, the global label gives a
        // consistent identity even though the same ancestral copy may now sit at
        // different sorted indices in different descendants.
        //
        // label → Vec<(parent_copy_nr, direction)> — one entry per descendant
        let mut label_info: HashMap<usize, Vec<(usize, &'static str)>> = HashMap::new();
        let mut n_contributing = 0usize;

        for &desc_id in &known_desc {
            let desc_hap = match hap_lookup(base_haps, sim_haps, desc_id) { Some(h) => h, None => continue };
            let desc_alleles: Vec<_> = desc_hap.get_alleles_by_marker_name(&marker.name);
            if desc_alleles.len() != n_copies { continue; }

            let matrix = SimpleDifferenceMatrix::new(parent_alleles.clone(), desc_alleles);
            let deltas = matrix.calculate_mutations();

            let desc_labels = labeling.get_labels(desc_id, &marker.name);

            for delta in &deltas {
                let dir: &'static str = if delta.step > 0 { "up" } else if delta.step < 0 { "down" } else { "none" };
                // Label = global label of the descendant's copy, or fallback to child_copy index
                let label = desc_labels
                    .and_then(|lv| lv.get(delta.child_copy as usize).copied())
                    .unwrap_or(delta.child_copy as usize);
                label_info.entry(label).or_default().push((delta.parent_copy as usize, dir));
            }
            n_contributing += 1;
        }

        if n_contributing == 0 { continue; }

        for (label, infos) in &label_info {
            // Require every contributing descendant to have provided a direction for this label
            if infos.len() < n_contributing { continue; }

            let first_dir = infos[0].1;
            if first_dir == "none" { continue; }
            let all_agree = infos.iter().all(|(_, d)| *d == first_dir);
            if !all_agree { continue; }

            let direction = match first_dir {
                "up"   => BiasDirection::Up,
                "down" => BiasDirection::Down,
                _      => continue,
            };

            // Determine which parent copy to bias via majority vote across descendants
            let mut copy_counts: HashMap<usize, usize> = HashMap::new();
            for (pc, _) in infos { *copy_counts.entry(*pc).or_default() += 1; }
            let parent_copy_nr = copy_counts.into_iter()
                .max_by_key(|(_, c)| *c)
                .map(|(k, _)| k)
                .unwrap_or(*label % n_copies);

            let target_mass = match bias_value {
                Some(bv) => bv,
                None => (auto_bias_strength / (1.0 + distance_to_mrca as f64))
                    .max(auto_bias_min)
                    .min(auto_bias_max),
            };
            biases.push(Bias {
                marker: marker.clone(),
                copy_nr: parent_copy_nr as u32,
                direction,
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

    // Pre-compute globally consistent copy labels for all known individuals
    let copy_labeling = CopyLabeling::compute(pedigree, marker_set);

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

    // Debug: write a diagnostic if any known→known edge already has zero probability.
    // This means ALL iterations will return P=0 regardless of how unknowns are sampled.
    let debug_max = params.debug_zero_prob_samples.unwrap_or(0) as usize;
    if debug_max > 0 && !log_const_edge_prob.is_finite() {
        use std::fmt::Write as FmtWrite;
        let mut out = String::new();
        let _ = writeln!(out, "MatchY — Zero-Probability Debug");
        let _ = writeln!(out, "\nCause: a known→known edge has P=0 (requires >2-step mutation).");
        let _ = writeln!(out, "This makes ALL simulation iterations return P=0.\n");
        let _ = writeln!(out, "Known→known edge probabilities:");
        for r in pedigree.relationships.iter().filter(|r| {
            !unknown_id_set.contains(r.parent_id.as_str())
                && !unknown_id_set.contains(r.child_id.as_str())
        }) {
            if let (Some(ph), Some(ch)) = (
                base_haplotypes.get(&r.parent_id),
                base_haplotypes.get(&r.child_id),
            ) {
                let p = get_edge_probability(ph, ch, marker_set, params.two_step_mutation_fraction);
                let flag = if p == 0.0 { "  ← ZERO" } else { "" };
                let _ = writeln!(out, "  {} → {}:  p={:.4e}{}", r.parent_id, r.child_id, p, flag);
                if p == 0.0 {
                    for m in &marker_set.markers {
                        let pa = ph.get_alleles_by_marker_name(&m.name);
                        let ca = ch.get_alleles_by_marker_name(&m.name);
                        if pa.is_empty() || ca.is_empty() || pa.len() != ca.len() { continue; }
                        let mp = calculate_mutation_probability(
                            &pa, &ca, m.single_copy_mutation_rate(), params.two_step_mutation_fraction,
                        );
                        if mp == 0.0 {
                            let pstr: Vec<String> = pa.iter().map(|a| a.display()).collect();
                            let cstr: Vec<String> = ca.iter().map(|a| a.display()).collect();
                            let _ = writeln!(
                                out,
                                "    {:20} parent=[{}]  child=[{}]  ← ZERO",
                                m.name, pstr.join(", "), cstr.join(", ")
                            );
                        }
                    }
                }
            }
        }
        let debug_dir = params.debug_zero_prob_path.as_deref().unwrap_or(&params.results_path);
        let debug_path = debug_dir.join("debug_zero_prob.txt");
        if let Some(parent) = debug_path.parent() {
            if !parent.as_os_str().is_empty() {
                let _ = std::fs::create_dir_all(parent);
            }
        }
        let _ = std::fs::write(&debug_path, out);
        tracing::info!("Wrote zero-prob const-edge diagnostic to {}", debug_path.display());
    }

    // Debug: collector for per-iteration zero-probability samples (variable edges).
    let debug_collector: Option<std::sync::Arc<std::sync::Mutex<Vec<ZeroProbSample>>>> =
        if debug_max > 0 {
            Some(std::sync::Arc::new(std::sync::Mutex::new(Vec::with_capacity(debug_max))))
        } else {
            None
        };

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
                    // Log-space accumulation avoids underflow for large pedigrees.
                    // log_iw = sum(ln(u_i) - ln(w_i)) over all allele copies of all unknowns.
                    let mut log_iw_total = 0.0f64;

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
                            pedigree, uid, marker_set, &base_haplotypes, &sim,
                            params.bias, params.auto_bias_strength, params.auto_bias_min, params.auto_bias_max,
                            None, &copy_labeling,
                        );

                        // (C) Use precomputed neutral step tables
                        let (new_hap, w, u) = mutate_haplotype_precomputed(
                            &parent_hap, marker_set, &neutral_tables, &biases, &mut local_rng,
                        );

                        sim.insert(uid.clone(), new_hap);
                        log_iw_total += u.ln() - w.ln();
                    }

                    // (D+5) Known→known edges: -inf means at least one probability is zero.
                    if !log_const_edge_prob.is_finite() {
                        return (0.0f64, log_iw_total);
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
                        if let Some(ref collector) = debug_collector {
                            let mut guard = collector.lock().unwrap();
                            if guard.len() < debug_max {
                                let sample = collect_zero_prob_sample(
                                    &sim, &base_haplotypes, &ordered_unknown, &child_of,
                                    &var_edge_rels, marker_set, params.two_step_mutation_fraction,
                                );
                                guard.push(sample);
                            }
                        }
                        return (0.0f64, log_iw_total);
                    }

                    let ped_prob = (log_const_edge_prob + log_var_edge).exp();
                    let probability = if ped_prob.is_finite() { ped_prob } else { 0.0 };
                    // Return (probability, log_iw_total) — linear conversion happens after
                    // batch-max normalization below.
                    (probability, log_iw_total)
                })
                .collect()
        })
        .collect();

    // Normalize importance weights by the batch maximum before converting from log-space.
    // The self-normalizing IS ratio sum(p·iw)/sum(iw) is invariant to multiplying all
    // iw by a constant, so subtracting max_log_iw cancels in the ratio while preventing
    // exp() from underflowing to 0 regardless of pedigree size.
    let max_log_iw: f64 = all_samples.iter().flatten()
        .map(|&(_, l)| l)
        .fold(f64::NEG_INFINITY, f64::max);

    for chunk in all_samples {
        for (probability, log_iw) in chunk {
            let importance_weight = if max_log_iw.is_finite() {
                (log_iw - max_log_iw).exp()
            } else {
                0.0
            };
            result.accumulate(probability, importance_weight);
        }
    }

    // Write debug samples collected during this batch (only if we captured anything new)
    if let Some(collector) = debug_collector {
        let samples = collector.lock().unwrap();
        if !samples.is_empty() {
            let debug_dir = params.debug_zero_prob_path.as_deref().unwrap_or(&params.results_path);
            let debug_path = debug_dir.join("debug_zero_prob.txt");
            match write_zero_prob_debug(&samples, &debug_path) {
                Ok(()) => tracing::info!(
                    "Wrote {} zero-probability debug sample(s) to {}",
                    samples.len(),
                    debug_path.display()
                ),
                Err(e) => tracing::warn!("Failed to write zero-prob debug file: {}", e),
            }
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
    avg_pedigree_probability: f64,
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

    // Pre-compute globally consistent copy labels for all known individuals
    let copy_labeling_match = CopyLabeling::compute(pedigree, marker_set);

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
    let avg_pp = avg_pedigree_probability;

    // One RNG per worker thread — avoids creating batch_length separate RNG objects
    let n_chunks = rayon::current_num_threads().max(1);
    let chunk_size = (batch_length as usize + n_chunks - 1) / n_chunks;
    let chunk_seeds: Vec<u64> = (0..n_chunks).map(|_| rng.gen()).collect();

    // Each item: (probability, log_iw, match_count_opt, matching_individual_ids).
    // Weighted deltas are computed after batch-max normalization, not inside the parallel closure.
    type IterSample = (f64, f64, Option<u32>, Vec<String>);

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

                    // Log-space accumulation avoids underflow for large pedigrees.
                    let mut log_iw_total = 0.0f64;
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
                            pedigree, uid, marker_set, &base_haplotypes, &sim,
                            params.bias, params.auto_bias_strength, params.auto_bias_min, params.auto_bias_max,
                            Some(fixed_id.as_str()), &copy_labeling_match,
                        );
                        // (C) Use precomputed neutral step tables
                        let (new_hap, w, u) = mutate_haplotype_precomputed(
                            &parent_hap, marker_set, &neutral_tables_match, &biases, &mut local_rng,
                        );
                        sim.insert(uid.clone(), new_hap);
                        log_iw_total += u.ln() - w.ln();
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

                    let probability = if probability.is_finite() { probability } else { 0.0 };

                    let match_count = if !is_outside && !non_excl_matches.is_empty() {
                        Some(non_excl_matches.len() as u32)
                    } else {
                        None
                    };

                    // Return raw log_iw; linear conversion + weighted deltas happen after
                    // batch-max normalization in the accumulation pass below.
                    (probability, log_iw_total, match_count, non_excl_matches)
                })
                .collect()
        })
        .collect();

    // Normalize by batch maximum before converting from log-space (same as pedigree batch).
    let max_log_iw: f64 = all_samples.iter().flatten()
        .map(|&(_, l, _, _)| l)
        .fold(f64::NEG_INFINITY, f64::max);

    for chunk in all_samples {
        for (prob, log_iw, match_count, matching_inds) in chunk {
            let iw = if max_log_iw.is_finite() {
                (log_iw - max_log_iw).exp()
            } else {
                0.0
            };
            let weighted = prob * iw;
            if let Some(k) = match_count {
                *result.match_accumulators.entry(k).or_default() += weighted;
            }
            for id in matching_inds {
                *result.per_individual.entry(id).or_default() += weighted;
            }
            result.accumulate(prob, iw);
        }
    }

    Ok(result)
}
