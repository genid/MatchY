/// Three-model ensemble convergence loop.
///
/// Runs 3 independent Monte Carlo models. After each batch, checks whether
/// all 3 model means are within `convergence_criterion` of their grand mean.
/// On divergence, retries up to MAX_TRIALS times.
///
/// Parallelism: the 3 models within each trial run concurrently via Rayon.
/// Each model gets its own seeded RNG so results are deterministic.
use crate::{MatchyError, Pedigree, SimulationParameters};
use crate::simulation::bias::AdaptiveBiasSchedule;
use crate::simulation::importance::{
    BatchResult,
    simulate_pedigree_probability_batch,
    simulate_matching_haplotypes_batch,
};
use crate::simulation::ProgressEvent;
use crate::{Haplotype, MarkerSet};
use crate::Result;
use rand::SeedableRng;
use rand::rngs::StdRng;
use rayon::prelude::*;
use rust_decimal::Decimal;
use std::collections::HashMap;

const NUM_MODELS: usize = 3;
const MAX_TRIALS: u32 = 10;

// ---------------------------------------------------------------------------
// check_convergence
// ---------------------------------------------------------------------------

/// Check if all 3 model means are within `criterion` of their grand mean.
pub fn check_convergence(means: &[Decimal; NUM_MODELS], criterion: f64) -> bool {
    let grand_mean: Decimal =
        means.iter().sum::<Decimal>() / Decimal::from(NUM_MODELS as u32);
    if grand_mean == Decimal::ZERO {
        return false;
    }
    let criterion_dec = Decimal::try_from(criterion).unwrap_or(Decimal::new(2, 2));
    means.iter().all(|m| {
        let diff = (*m - grand_mean).abs();
        diff / grand_mean <= criterion_dec
    })
}

// ---------------------------------------------------------------------------
// EnsembleTrial
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct EnsembleTrial {
    pub trial_nr: u32,
    pub model_results: [BatchResult; NUM_MODELS],
    pub converged: bool,
    /// Grand mean after convergence
    pub grand_mean: Option<Decimal>,
}

impl EnsembleTrial {
    pub fn new(trial_nr: u32) -> Self {
        Self {
            trial_nr,
            model_results: [
                BatchResult::new(),
                BatchResult::new(),
                BatchResult::new(),
            ],
            converged: false,
            grand_mean: None,
        }
    }

    pub fn running_means(&self) -> [Decimal; NUM_MODELS] {
        std::array::from_fn(|i| {
            self.model_results[i]
                .running_mean()
                .unwrap_or(Decimal::ZERO)
        })
    }

    pub fn grand_mean(&self) -> Decimal {
        let means = self.running_means();
        means.iter().sum::<Decimal>() / Decimal::from(NUM_MODELS as u32)
    }

    pub fn check_convergence(&self, criterion: f64) -> bool {
        let means = self.running_means();
        check_convergence(&means, criterion)
    }
}

// ---------------------------------------------------------------------------
// run_ensemble  (pedigree probability step — Step 1)
// ---------------------------------------------------------------------------

/// Run the Step 1 ensemble: average pedigree probability.
///
/// Returns (grand_mean, per_model_means, trial_count).
pub fn run_ensemble_pedigree_probability(
    pedigree: &Pedigree,
    root_id: &str,
    marker_set: &MarkerSet,
    params: &SimulationParameters,
    progress_tx: Option<&std::sync::mpsc::Sender<ProgressEvent>>,
    adaptive_schedule: Option<&mut AdaptiveBiasSchedule>,
) -> Result<EnsembleTrial> {
    use crate::simulation::SimulationStage;

    let mut trial = EnsembleTrial::new(1);
    let mut tightening = 0u32;
    let mut current_criterion = params.convergence_criterion;
    let mut adaptive = adaptive_schedule;

    for trial_nr in 1..=MAX_TRIALS {
        trial = EnsembleTrial::new(trial_nr);

        // Determine per-model bias values (adaptive or fixed)
        let model_biases: [Option<f64>; NUM_MODELS] = if params.adaptive_bias {
            if let Some(ref sched) = adaptive {
                sched.model_biases.map(Some)
            } else {
                [Some(0.10); NUM_MODELS]
            }
        } else {
            [params.bias; NUM_MODELS]
        };

        // Run the 3 models in parallel. Each has its own seeded RNG and cloned
        // params so there is no shared mutable state.
        let model_batches: Vec<Result<BatchResult>> = (0..NUM_MODELS)
            .into_par_iter()
            .map(|model| {
                let mut model_params = params.clone();
                model_params.bias = model_biases[model];
                let seed = (trial_nr as u64) * 100 + model as u64;
                let mut rng = StdRng::seed_from_u64(seed);
                simulate_pedigree_probability_batch(
                    pedigree,
                    root_id,
                    marker_set,
                    &model_params,
                    &mut rng,
                    params.batch_length,
                )
            })
            .collect();

        // Propagate first error (if any), then emit progress and store results.
        for (model, batch_result) in model_batches.into_iter().enumerate() {
            let batch = batch_result?;
            if let Some(tx) = progress_tx {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: model as u8,
                    iteration: batch.iterations,
                    current_mean: batch
                        .running_mean()
                        .map(|m| m.to_string())
                        .unwrap_or_else(|| "0".into()),
                    stage: SimulationStage::PedigreeProbability,
                    converged: false,
                });
            }
            trial.model_results[model] = batch;
        }

        let means = trial.running_means();
        let grand = trial.grand_mean();

        // Check if grand mean > 1 (invalid — tighten criterion)
        if grand > Decimal::ONE {
            tracing::warn!(
                "Trial {}: grand mean {} > 1. Tightening criterion to {}",
                trial_nr,
                grand,
                current_criterion / 1.5
            );
            if tightening >= 2 {
                return Err(MatchyError::SimulationFailed(
                    "Average pedigree probability exceeds 1 after max tightening".into(),
                ));
            }
            current_criterion /= 1.5;
            tightening += 1;
            continue;
        }

        if check_convergence(&means, current_criterion) {
            trial.converged = true;
            trial.grand_mean = Some(grand);
            tracing::info!(
                "Pedigree probability converged after {} trials: {}",
                trial_nr,
                grand
            );

            if let Some(ref mut sched) = adaptive {
                let estimates = means.map(|m| f64::try_from(m).unwrap_or(0.0));
                sched.update(estimates);
            }

            // Emit final convergence event
            if let Some(tx) = progress_tx {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: 0,
                    iteration: trial.model_results[0].iterations,
                    current_mean: grand.to_string(),
                    stage: SimulationStage::PedigreeProbability,
                    converged: true,
                });
            }

            return Ok(trial);
        }

        tracing::info!(
            "Trial {} did not converge. Means: {:?}. Grand: {}",
            trial_nr,
            means,
            grand
        );
    }

    Err(MatchyError::ConvergenceFailed(format!(
        "Pedigree probability did not converge after {} trials",
        MAX_TRIALS
    )))
}

// ---------------------------------------------------------------------------
// run_ensemble_matching_haplotypes  (Step 2)
// ---------------------------------------------------------------------------

/// Run the Step 2 ensemble: inside-pedigree match probabilities.
///
/// Returns EnsembleTrial with match_accumulators and per_individual filled.
pub fn run_ensemble_matching_haplotypes(
    pedigree: &Pedigree,
    root_id: &str,
    suspect_haplotype: &Haplotype,
    avg_pedigree_probability: Decimal,
    marker_set: &MarkerSet,
    params: &SimulationParameters,
    is_outside: bool,
    progress_tx: Option<&std::sync::mpsc::Sender<ProgressEvent>>,
    adaptive_schedule: Option<&mut AdaptiveBiasSchedule>,
) -> Result<EnsembleTrial> {
    use crate::simulation::SimulationStage;

    let stage = if is_outside {
        SimulationStage::OutsideMatchProbability
    } else {
        SimulationStage::InsideMatchProbability
    };

    let mut current_criterion = params.convergence_criterion;
    let mut tightening = 0u32;
    let mut adaptive = adaptive_schedule;

    for trial_nr in 1..=MAX_TRIALS {
        let mut trial = EnsembleTrial::new(trial_nr);

        let model_biases: [Option<f64>; NUM_MODELS] = if params.adaptive_bias {
            if let Some(ref sched) = adaptive {
                sched.model_biases.map(Some)
            } else {
                [Some(0.10); NUM_MODELS]
            }
        } else {
            [params.bias; NUM_MODELS]
        };

        // Run the 3 models in parallel.
        let model_batches: Vec<Result<BatchResult>> = (0..NUM_MODELS)
            .into_par_iter()
            .map(|model| {
                let mut model_params = params.clone();
                model_params.bias = model_biases[model];
                let seed = (trial_nr as u64 + 1000) * 100 + model as u64;
                let mut rng = StdRng::seed_from_u64(seed);
                simulate_matching_haplotypes_batch(
                    pedigree,
                    root_id,
                    suspect_haplotype,
                    avg_pedigree_probability,
                    marker_set,
                    &model_params,
                    is_outside,
                    &mut rng,
                    params.batch_length,
                )
            })
            .collect();

        for (model, batch_result) in model_batches.into_iter().enumerate() {
            let batch = batch_result?;
            if let Some(tx) = progress_tx {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: model as u8,
                    iteration: batch.iterations,
                    current_mean: batch
                        .running_mean()
                        .map(|m| m.to_string())
                        .unwrap_or_else(|| "0".into()),
                    stage,
                    converged: false,
                });
            }
            trial.model_results[model] = batch;
        }

        let means = trial.running_means();
        let grand = trial.grand_mean();

        if check_convergence(&means, current_criterion) {
            trial.converged = true;
            trial.grand_mean = Some(grand);

            if let Some(ref mut sched) = adaptive {
                let estimates = means.map(|m| f64::try_from(m).unwrap_or(0.0));
                sched.update(estimates);
            }

            if let Some(tx) = progress_tx {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: 0,
                    iteration: trial.model_results[0].iterations,
                    current_mean: grand.to_string(),
                    stage,
                    converged: true,
                });
            }

            return Ok(trial);
        }

        tracing::info!(
            "Match probability trial {} did not converge. Grand: {}",
            trial_nr, grand
        );
    }

    Err(MatchyError::ConvergenceFailed(format!(
        "Match probability did not converge after {} trials",
        MAX_TRIALS
    )))
}

// ---------------------------------------------------------------------------
// aggregate_ensemble_result
// ---------------------------------------------------------------------------

/// Aggregate per-individual and per-match-count probabilities from 3 models.
pub fn aggregate_per_individual(models: &[BatchResult; NUM_MODELS]) -> HashMap<String, Decimal> {
    let mut all_ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    for m in models {
        all_ids.extend(m.per_individual.keys().cloned());
    }

    let mut result = HashMap::new();
    for id in all_ids {
        let sum: Decimal = models
            .iter()
            .map(|m| {
                let w_sum = m.per_individual.get(&id).copied().unwrap_or(Decimal::ZERO);
                let denom = m.weight_sum;
                if denom == Decimal::ZERO { Decimal::ZERO } else { w_sum / denom }
            })
            .sum();
        result.insert(id, sum / Decimal::from(NUM_MODELS as u32));
    }
    result
}

/// Aggregate per-match-count probabilities from 3 models.
pub fn aggregate_match_counts(models: &[BatchResult; NUM_MODELS]) -> HashMap<u32, Decimal> {
    let mut all_counts: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for m in models {
        all_counts.extend(m.match_accumulators.keys().copied());
    }

    let mut result = HashMap::new();
    for count in all_counts {
        let sum: Decimal = models
            .iter()
            .map(|m| {
                let w_sum = m.match_accumulators.get(&count).copied().unwrap_or(Decimal::ZERO);
                let denom = m.weight_sum;
                if denom == Decimal::ZERO { Decimal::ZERO } else { w_sum / denom }
            })
            .sum();
        result.insert(count, sum / Decimal::from(NUM_MODELS as u32));
    }
    result
}
