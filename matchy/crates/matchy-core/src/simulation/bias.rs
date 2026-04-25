/// Importance sampling bias calculation.
///
/// For each unknown individual, compute the Bias structs that steer mutation
/// sampling toward/away from matching the suspect's haplotype.
///
/// The bias formula (default, when params.bias is None):
///   target_mass = min(max(0.1, 0.8 / (1 + distance_to_mrca)), 0.4)
///
/// In trace mode, no bias is applied (uniform picking).
use crate::Bias;
use rustc_hash::FxHashMap;

pub type BiasCache = FxHashMap<u64, Vec<Bias>>;

/// Compute the auto-bias target mass for an individual at a given distance
/// from the MRCA of its known descendants.
///
/// Formula: min(max(0.1, 0.8 / (1 + d)), 0.4) — mirrors Python models.py:914.
pub fn default_bias_target_mass(distance_to_mrca: u32) -> f64 {
    (0.8 / (1.0 + distance_to_mrca as f64))
        .max(0.1)
        .min(0.4)
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
        // Rank models: best = highest mean (mirrors Python: sorted by -x[1]).
        // Best model gets the lowest bias (0.05); worst gets highest (0.25).
        let mut ranked: Vec<(usize, f64)> = estimates.iter().copied().enumerate().collect();
        ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        self.model_biases[ranked[0].0] = 0.05; // best
        self.model_biases[ranked[1].0] = 0.15; // middle
        self.model_biases[ranked[2].0] = 0.25; // worst
    }
}
