/// Importance sampling bias calculation.
///
/// For each unknown individual, compute the Bias structs that steer mutation
/// sampling toward/away from matching the suspect's haplotype.
///
/// The bias formula (default, when params.bias is None):
///   target_mass = min(max(0.1, 0.8 / (1 + distance_to_mrca)), 0.4)
///
/// In trace mode, no bias is applied (uniform picking).
use crate::{Bias, BiasDirection, Haplotype, Individual, Marker, MarkerSet, Pedigree};
use crate::hungarian::SimpleDifferenceMatrix;
use rustc_hash::FxHashMap;

/// Cache key for get_biases: hash of the pedigree + suspect haplotype state.
/// In Python this is handled by lru_cache keyed on (sorted individual haplotypes).
pub type BiasCache = FxHashMap<u64, Vec<Bias>>;

/// Compute the auto-bias target mass for an individual at a given distance
/// from the MRCA.
pub fn default_bias_target_mass(distance_to_mrca: u32) -> f64 {
    (0.8 / (1.0 + distance_to_mrca as f64))
        .max(0.1)
        .min(0.4)
}

/// Compute Bias structs for a given individual relative to the suspect haplotype.
///
/// For each marker, uses the Hungarian algorithm to find the optimal parent→child
/// allele assignment and derives a Bias toward the direction that moves the
/// individual's alleles toward the suspect's alleles.
pub fn compute_biases_for_individual(
    individual: &Individual,
    suspect_haplotype: &Haplotype,
    marker_set: &MarkerSet,
    distance_to_mrca: u32,
    fixed_bias: Option<f64>,
) -> Vec<Bias> {
    let target_mass = fixed_bias.unwrap_or_else(|| default_bias_target_mass(distance_to_mrca));

    let mut biases = Vec::new();
    for marker in &marker_set.markers {
        let individual_alleles = individual.get_alleles_by_marker_name(&marker.name);
        let suspect_alleles = suspect_haplotype.get_alleles_by_marker_name(&marker.name);

        if individual_alleles.is_empty() || suspect_alleles.is_empty() {
            continue;
        }
        if individual_alleles.len() != suspect_alleles.len() {
            continue;
        }

        // Cast to owned refs for the matrix
        let parent_alleles: Vec<&crate::Allele> = individual_alleles;
        let child_alleles: Vec<&crate::Allele> = suspect_alleles;

        let matrix = SimpleDifferenceMatrix::new(parent_alleles, child_alleles);
        let deltas = matrix.calculate_mutations();

        for (copy_nr, delta) in deltas.iter().enumerate() {
            if delta.step == 0 {
                continue; // Already matching — no bias needed
            }
            let direction = if delta.step > 0 {
                BiasDirection::Up
            } else {
                BiasDirection::Down
            };
            biases.push(Bias {
                marker: marker.clone(),
                copy_nr: copy_nr as u32,
                direction,
                target_mass,
            });
        }
    }
    biases
}

/// Adaptive bias schedule for multi-trial convergence.
///
/// After each trial, models are ranked by their last estimated mean probability.
/// Best model gets the lowest bias (0.05), middle gets 0.15, worst gets 0.25.
/// Trial 1 always uses 0.10 for all models.
#[derive(Debug, Clone)]
pub struct AdaptiveBiasSchedule {
    pub trial: u32,
    /// Per-model bias values for the next trial
    pub model_biases: [f64; 3],
}

impl Default for AdaptiveBiasSchedule {
    fn default() -> Self {
        Self {
            trial: 0,
            model_biases: [0.10, 0.10, 0.10],
        }
    }
}

impl AdaptiveBiasSchedule {
    /// Update the schedule based on the last estimates from each model.
    /// `estimates[i]` is the last running mean for model i.
    pub fn update(&mut self, estimates: [f64; 3]) {
        self.trial += 1;
        // Rank models: best estimate = lowest absolute value (closest to convergence target)
        let mut ranked: Vec<(usize, f64)> = estimates.iter().copied().enumerate().collect();
        ranked.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        self.model_biases[ranked[0].0] = 0.05; // best
        self.model_biases[ranked[1].0] = 0.15; // middle
        self.model_biases[ranked[2].0] = 0.25; // worst
    }
}
