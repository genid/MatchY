/// Three-model ensemble convergence loop.
///
/// Runs 3 independent Monte Carlo models. After each batch, checks whether
/// all 3 model means are within `convergence_criterion` of their grand mean.
/// Keeps running until convergence is reached; there is no trial limit.
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
use std::sync::atomic::{AtomicBool, Ordering};

const NUM_MODELS: usize = 3;

// ---------------------------------------------------------------------------
// check_convergence
// ---------------------------------------------------------------------------

/// Check convergence by requiring every per-iteration running mean across all 3
/// models to be within `criterion` of the final grand mean.
///
/// Mirrors Python simulation.py:321–326 exactly:
///   total_mean = mean(model_probabilities[m][-1] for m in range(3))
///   if all(all(abs(prob - total_mean) / total_mean < threshold
///              for prob in model_probabilities[m]) for m in range(3)):
///
/// `BatchResult::running_means` plays the role of Python's `model_probabilities[m]`:
/// it is a list of per-iteration running means, optionally prefixed with the last
/// mean from the previous trial (the `model_probabilities[m][-1:]` carryover).
pub fn check_convergence(results: &[BatchResult; NUM_MODELS], criterion: f64) -> bool {
    let final_means: Vec<Decimal> = results
        .iter()
        .filter_map(|r| r.running_means.last().copied())
        .collect();
    if final_means.len() < NUM_MODELS {
        return false;
    }
    let grand_mean: Decimal =
        final_means.iter().sum::<Decimal>() / Decimal::from(NUM_MODELS as u32);
    if grand_mean == Decimal::ZERO {
        return false;
    }
    let criterion_dec = Decimal::try_from(criterion).unwrap_or(Decimal::new(2, 2));
    results.iter().all(|r| {
        r.running_means.iter().all(|&prob| {
            let diff = (prob - grand_mean).abs();
            diff / grand_mean <= criterion_dec
        })
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

    /// Running means of each model (final weighted mean).
    pub fn running_means(&self) -> [Decimal; NUM_MODELS] {
        std::array::from_fn(|i| {
            self.model_results[i]
                .running_mean()
                .unwrap_or(Decimal::ZERO)
        })
    }

    /// Grand mean computed from the last entry of each model's running_means history
    /// (equivalent to the mean of the 3 models' final weighted means).
    pub fn grand_mean(&self) -> Decimal {
        let final_means: Vec<Decimal> = self.model_results.iter()
            .filter_map(|r| r.running_means.last().copied())
            .collect();
        if final_means.len() < NUM_MODELS {
            return Decimal::ZERO;
        }
        final_means.iter().sum::<Decimal>() / Decimal::from(NUM_MODELS as u32)
    }

    /// Extract the last running mean from each model for seeding the next trial.
    /// Mirrors Python: `model_probabilities = {m: model_probabilities[m][-1:] ...}`.
    pub fn last_means(&self) -> [Option<Decimal>; NUM_MODELS] {
        std::array::from_fn(|i| self.model_results[i].running_means.last().copied())
    }
}

// ---------------------------------------------------------------------------
// run_ensemble  (pedigree probability step — Step 1)
// ---------------------------------------------------------------------------

/// Run the Step 1 ensemble: average pedigree probability.
pub fn run_ensemble_pedigree_probability(
    pedigree: &Pedigree,
    root_id: &str,
    marker_set: &MarkerSet,
    params: &SimulationParameters,
    progress_tx: Option<&std::sync::mpsc::Sender<ProgressEvent>>,
    adaptive_schedule: Option<&mut AdaptiveBiasSchedule>,
    cancel_flag: Option<&AtomicBool>,
    stage: crate::simulation::SimulationStage,
) -> Result<EnsembleTrial> {
    let mut tightening = 0u32;
    let mut current_criterion = params.convergence_criterion;
    let mut adaptive = adaptive_schedule;
    let mut seeds: [Option<Decimal>; NUM_MODELS] = [None; NUM_MODELS];
    // Cumulative sums carried across trials — mirrors Python where weight_sums /
    // weighted_sums are never reset between trials.
    let mut carry_sums: [(f64, f64); NUM_MODELS] = [(0.0, 0.0); NUM_MODELS];

    for trial_nr in 1u32.. {
        if cancel_flag.map(|f| f.load(Ordering::Relaxed)).unwrap_or(false) {
            return Err(MatchyError::Cancelled);
        }

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

        let carry = carry_sums;
        let model_batches: Vec<Result<BatchResult>> = (0..NUM_MODELS)
            .into_par_iter()
            .map(|model| {
                let mut model_params = params.clone();
                model_params.bias = model_biases[model];
                let seed = params.seed.unwrap_or(0).wrapping_add((trial_nr as u64) * 100 + model as u64);
                let mut rng = StdRng::seed_from_u64(seed);
                simulate_pedigree_probability_batch(
                    pedigree,
                    root_id,
                    marker_set,
                    &model_params,
                    &mut rng,
                    params.batch_length,
                    Some(carry[model]),
                )
            })
            .collect();

        for (model, batch_result) in model_batches.into_iter().enumerate() {
            let mut batch = batch_result?;
            if let Some(seed) = seeds[model] {
                batch.running_means.insert(0, seed);
            }
            trial.model_results[model] = batch;
        }

        let grand = trial.grand_mean();

        if grand > Decimal::ONE {
            tracing::warn!(
                "Trial {}: grand mean {} > 1. Tightening criterion to {}",
                trial_nr,
                grand,
                current_criterion / 1.5
            );
            if let Some(tx) = progress_tx {
                for model in 0..NUM_MODELS {
                    let _ = tx.send(ProgressEvent {
                        trial: trial_nr,
                        model: model as u8,
                        iteration: trial.model_results[model].iterations,
                        current_mean: trial.model_results[model]
                            .running_mean()
                            .map(|m| m.to_string())
                            .unwrap_or_else(|| "0".into()),
                        stage,
                        converged: false,
                    });
                }
            }
            if tightening >= 2 {
                return Err(MatchyError::SimulationFailed(
                    "Average pedigree probability exceeds 1 after max tightening".into(),
                ));
            }
            seeds = trial.last_means();
            for m in 0..NUM_MODELS {
                carry_sums[m] = (trial.model_results[m].weighted_sum, trial.model_results[m].weight_sum);
            }
            current_criterion /= 1.5;
            tightening += 1;
            continue;
        }

        let converged_now = check_convergence(&trial.model_results, current_criterion);

        // Send progress events for all 3 models with the correct converged flag.
        // Doing this after the convergence check ensures all models carry the same
        // converged=true signal in the same batch — no extra event for model 0 only.
        if let Some(tx) = progress_tx {
            for model in 0..NUM_MODELS {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: model as u8,
                    iteration: trial.model_results[model].iterations,
                    current_mean: trial.model_results[model]
                        .running_mean()
                        .map(|m| m.to_string())
                        .unwrap_or_else(|| "0".into()),
                    stage,
                    converged: converged_now,
                });
            }
        }

        if converged_now {
            trial.converged = true;
            trial.grand_mean = Some(grand);
            tracing::info!(
                "Pedigree probability converged after {} trials: {}",
                trial_nr,
                grand
            );

            if let Some(ref mut sched) = adaptive {
                let estimates = trial.running_means().map(|m| f64::try_from(m).unwrap_or(0.0));
                sched.update(estimates);
            }

            return Ok(trial);
        }

        // Update adaptive schedule so next trial uses differentiated model biases.
        // Mirrors Python simulation.py:244-258 where per_model_bias is recomputed
        // at the START of each trial >= 2 based on the previous trial's last means.
        if let Some(ref mut sched) = adaptive {
            let estimates = trial.running_means().map(|m| f64::try_from(m).unwrap_or(0.0));
            sched.update(estimates);
        }

        seeds = trial.last_means();
        for m in 0..NUM_MODELS {
            carry_sums[m] = (trial.model_results[m].weighted_sum, trial.model_results[m].weight_sum);
        }
        tracing::info!(
            "Trial {} did not converge. Grand: {}. Adaptive biases: {:?}",
            trial_nr,
            grand,
            adaptive.as_ref().map(|s| s.model_biases)
        );
    }

    unreachable!("unbounded convergence loop exited without returning")
}

// ---------------------------------------------------------------------------
// run_ensemble_matching_haplotypes  (Step 2 / Step 3)
// ---------------------------------------------------------------------------

/// Run the Step 2 (inside) or Step 3 (outside) ensemble: match probabilities.
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
    cancel_flag: Option<&AtomicBool>,
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
    let mut seeds: [Option<Decimal>; NUM_MODELS] = [None; NUM_MODELS];
    // Full BatchResult carry across trials — all accumulators (weighted_sum, weight_sum,
    // match_accumulators, per_individual) are cumulative, matching Python where these
    // are never reset between trials.
    let mut carry: [Option<BatchResult>; NUM_MODELS] = [None, None, None];

    for trial_nr in 1u32.. {
        if cancel_flag.map(|f| f.load(Ordering::Relaxed)).unwrap_or(false) {
            return Err(MatchyError::Cancelled);
        }

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

        let carry_snapshot: [Option<BatchResult>; NUM_MODELS] = std::array::from_fn(|i| carry[i].clone());
        let model_batches: Vec<Result<BatchResult>> = (0..NUM_MODELS)
            .into_par_iter()
            .map(|model| {
                let mut model_params = params.clone();
                model_params.bias = model_biases[model];
                let seed = params.seed.unwrap_or(0).wrapping_add((trial_nr as u64 + 1000) * 100 + model as u64);
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
                    carry_snapshot[model].clone(),
                )
            })
            .collect();

        for (model, batch_result) in model_batches.into_iter().enumerate() {
            let mut batch = batch_result?;
            if let Some(seed) = seeds[model] {
                batch.running_means.insert(0, seed);
            }
            trial.model_results[model] = batch;
        }

        let grand = trial.grand_mean();

        if grand > Decimal::ONE {
            tracing::warn!(
                "Match trial {}: grand mean {} > 1. Tightening criterion.",
                trial_nr, grand
            );
            if let Some(tx) = progress_tx {
                for model in 0..NUM_MODELS {
                    let _ = tx.send(ProgressEvent {
                        trial: trial_nr,
                        model: model as u8,
                        iteration: trial.model_results[model].iterations,
                        current_mean: trial.model_results[model]
                            .running_mean()
                            .map(|m| m.to_string())
                            .unwrap_or_else(|| "0".into()),
                        stage,
                        converged: false,
                    });
                }
            }
            if tightening >= 2 {
                return Err(MatchyError::SimulationFailed(
                    "Match probability exceeds 1 after max tightening — pedigree likely invalid".into(),
                ));
            }
            seeds = trial.last_means();
            for m in 0..NUM_MODELS {
                carry[m] = Some(trial.model_results[m].clone());
            }
            current_criterion /= 1.5;
            tightening += 1;
            continue;
        }

        let converged_now = check_convergence(&trial.model_results, current_criterion);

        if let Some(tx) = progress_tx {
            for model in 0..NUM_MODELS {
                let _ = tx.send(ProgressEvent {
                    trial: trial_nr,
                    model: model as u8,
                    iteration: trial.model_results[model].iterations,
                    current_mean: trial.model_results[model]
                        .running_mean()
                        .map(|m| m.to_string())
                        .unwrap_or_else(|| "0".into()),
                    stage,
                    converged: converged_now,
                });
            }
        }

        if converged_now {
            trial.converged = true;
            trial.grand_mean = Some(grand);
            tracing::info!(
                "Match probability converged after {} trials: {}",
                trial_nr, grand
            );

            if let Some(ref mut sched) = adaptive {
                let estimates = trial.running_means().map(|m| f64::try_from(m).unwrap_or(0.0));
                sched.update(estimates);
            }

            return Ok(trial);
        }

        if let Some(ref mut sched) = adaptive {
            let estimates = trial.running_means().map(|m| f64::try_from(m).unwrap_or(0.0));
            sched.update(estimates);
        }

        seeds = trial.last_means();
        for m in 0..NUM_MODELS {
            carry[m] = Some(trial.model_results[m].clone());
        }
        tracing::info!(
            "Match probability trial {} did not converge. Grand: {}. Adaptive biases: {:?}",
            trial_nr, grand,
            adaptive.as_ref().map(|s| s.model_biases)
        );
    }

    unreachable!("unbounded convergence loop exited without returning")
}

// ---------------------------------------------------------------------------
// aggregate helpers
// ---------------------------------------------------------------------------

/// Aggregate per-individual probabilities from 3 models.
pub fn aggregate_per_individual(models: &[BatchResult; NUM_MODELS]) -> HashMap<String, Decimal> {
    let mut all_ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    for m in models {
        all_ids.extend(m.per_individual.keys().cloned());
    }

    let mut result = HashMap::new();
    for id in all_ids {
        let sum: f64 = models
            .iter()
            .map(|m| {
                let w_sum = m.per_individual.get(&id).copied().unwrap_or(0.0);
                let denom = m.weight_sum;
                if denom == 0.0 { 0.0 } else { w_sum / denom }
            })
            .sum();
        let mean = sum / NUM_MODELS as f64;
        result.insert(id, Decimal::try_from(mean).unwrap_or(Decimal::ZERO));
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
        let sum: f64 = models
            .iter()
            .map(|m| {
                let w_sum = m.match_accumulators.get(&count).copied().unwrap_or(0.0);
                let denom = m.weight_sum;
                if denom == 0.0 { 0.0 } else { w_sum / denom }
            })
            .sum();
        let mean = sum / NUM_MODELS as f64;
        result.insert(count, Decimal::try_from(mean).unwrap_or(Decimal::ZERO));
    }
    result
}
