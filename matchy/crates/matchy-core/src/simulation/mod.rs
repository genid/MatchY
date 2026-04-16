pub mod bias;
pub mod convergence;
pub mod importance;
pub mod mutation;
pub mod probability;

use crate::{
    HaplotypeClass, MatchProbabilities, Pedigree, SimulationParameters, SimulationResult,
};
use crate::graph::{bfs_layers, validate_dag};
use crate::simulation::bias::AdaptiveBiasSchedule;
use crate::simulation::convergence::{
    aggregate_match_counts, aggregate_per_individual,
    run_ensemble_matching_haplotypes, run_ensemble_pedigree_probability,
};
use crate::Result;
use rust_decimal::Decimal;
use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Top-level entry point
// ---------------------------------------------------------------------------

/// Top-level simulation entry point.
///
/// Runs the full Monte Carlo simulation pipeline:
/// 1. Validate and prepare pedigree
/// 2. Run Step 1: average pedigree probability (P(hv))
/// 3. Run Step 2: inside-pedigree match probabilities (P(m|hv))
/// 4. Run Step 3: outside-pedigree match probability (P(m_outside))
/// 5. Return SimulationResult
pub fn run_simulation(
    mut pedigree: Pedigree,
    marker_set: &crate::MarkerSet,
    params: SimulationParameters,
    progress_tx: Option<std::sync::mpsc::Sender<ProgressEvent>>,
) -> Result<SimulationResult> {
    validate_dag(&pedigree)?;

    // Determine root (single root expected after pedigree preparation)
    let roots = pedigree.roots();
    let root_id = roots
        .first()
        .ok_or_else(|| crate::MatchyError::InvalidPedigree("No root individual".into()))?
        .to_string();

    // Find the suspect individual
    let suspect_haplotype = if let Some(ref suspect_name) = params.suspect {
        let ind = pedigree
            .get_individual_by_name(suspect_name)
            .ok_or_else(|| {
                crate::MatchyError::IndividualNotFound(suspect_name.clone())
            })?;
        Some(ind.haplotype.clone())
    } else {
        None
    };

    // Check if only root is known (shortcut)
    let layers = bfs_layers(&pedigree, &root_id)?;
    let bfs_order: Vec<String> = layers.into_iter().flatten().collect();
    let known_count = bfs_order
        .iter()
        .filter(|id| {
            pedigree
                .get_individual_by_id(id)
                .map(|i| i.haplotype_class != HaplotypeClass::Unknown)
                .unwrap_or(false)
        })
        .count();

    let mut adaptive = if params.adaptive_bias {
        Some(AdaptiveBiasSchedule::default())
    } else {
        None
    };

    // Step 1: Average pedigree probability
    let prog_ref = progress_tx.as_ref();
    let pedigree_prob_trial = if known_count <= 1 {
        tracing::info!("Only one known individual — pedigree probability is 1");
        // Return a trivial result
        let mut trial = convergence::EnsembleTrial::new(1);
        trial.converged = true;
        trial.grand_mean = Some(Decimal::ONE);
        trial
    } else {
        run_ensemble_pedigree_probability(
            &pedigree,
            &root_id,
            marker_set,
            &params,
            prog_ref,
            adaptive.as_mut(),
        )?
    };

    let avg_pedigree_probability = pedigree_prob_trial.grand_mean.unwrap_or(Decimal::ZERO);

    // Steps 2 and 3 only run if there is a suspect haplotype
    let inside_match_probabilities = if let Some(ref suspect_hap) = suspect_haplotype {
        if params.skip_inside {
            None
        } else {
            let trial = run_ensemble_matching_haplotypes(
                &pedigree,
                &root_id,
                suspect_hap,
                avg_pedigree_probability,
                marker_set,
                &params,
                false, // inside
                prog_ref,
                adaptive.as_mut(),
            )?;
            let probs = aggregate_match_counts(&trial.model_results);
            Some(MatchProbabilities {
                probabilities: probs,
                average_pedigree_probability: avg_pedigree_probability,
            })
        }
    } else {
        None
    };

    let outside_match_probability = if let Some(ref suspect_hap) = suspect_haplotype {
        if params.skip_outside {
            None
        } else {
            let trial = run_ensemble_matching_haplotypes(
                &pedigree,
                &root_id,
                suspect_hap,
                avg_pedigree_probability,
                marker_set,
                &params,
                true, // outside
                prog_ref,
                adaptive.as_mut(),
            )?;
            trial.grand_mean
        }
    } else {
        None
    };

    // Per-individual probabilities for trace mode
    let per_individual_probabilities = if params.trace_mode {
        if let Some(ref _trial) = inside_match_probabilities {
            // In trace mode we aggregate per-individual from the inside trial
            // Re-run won't be needed; it's populated from the inside run above
            None // TODO: plumb per_individual out from run_ensemble_matching_haplotypes
        } else {
            None
        }
    } else {
        None
    };

    Ok(SimulationResult {
        parameters: params,
        inside_match_probabilities,
        outside_match_probability,
        per_individual_probabilities,
        trials: pedigree_prob_trial.trial_nr,
        converged: pedigree_prob_trial.converged,
    })
}

// ---------------------------------------------------------------------------
// ProgressEvent (used by Tauri commands for live streaming)
// ---------------------------------------------------------------------------

/// Progress event emitted during simulation via mpsc channel.
#[derive(Debug, Clone, serde::Serialize)]
pub struct ProgressEvent {
    pub trial: u32,
    pub model: u8,
    pub iteration: u64,
    /// Running mean as string to avoid floating-point serialisation issues
    pub current_mean: String,
    pub stage: SimulationStage,
    pub converged: bool,
}

#[derive(Debug, Clone, Copy, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SimulationStage {
    PedigreeProbability,
    InsideMatchProbability,
    OutsideMatchProbability,
}
