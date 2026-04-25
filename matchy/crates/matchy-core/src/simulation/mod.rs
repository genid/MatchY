pub mod bias;
pub mod convergence;
pub mod importance;
pub mod mutation;
pub mod probability;

use crate::{
    HaplotypeClass, MatchProbabilities, Pedigree, SimulationParameters, SimulationResult,
    StageStats,
};
use crate::graph::{bfs_layers, validate_dag};
use crate::simulation::bias::AdaptiveBiasSchedule;
use crate::simulation::convergence::{
    aggregate_match_counts,
    aggregate_per_individual,
    run_ensemble_matching_haplotypes, run_ensemble_pedigree_probability,
    EnsembleTrial,
};
use crate::Result;
use rayon;
use rust_decimal::Decimal;
use std::sync::{Arc, atomic::{AtomicBool, Ordering}};
use std::time::Instant;

fn make_stage_stats(trial: &EnsembleTrial, batch_length: u64, elapsed: std::time::Duration) -> StageStats {
    let iters = trial.trial_nr as u64 * batch_length;
    StageStats {
        iterations_per_model: vec![iters; 3],
        model_probabilities: trial.model_results.iter()
            .filter_map(|br| br.running_means.last())
            .map(|d| format!("{:.4E}", f64::try_from(*d).unwrap_or(0.0)))
            .collect(),
        runtime_secs: elapsed.as_secs_f64(),
    }
}

// ---------------------------------------------------------------------------
// Top-level entry point
// ---------------------------------------------------------------------------

/// Top-level simulation entry point.
///
/// Configures a Rayon thread pool sized to `params.number_of_threads` and
/// delegates to `run_simulation_impl`. All nested `par_iter` calls (the 3-model
/// ensemble in `convergence.rs` and the per-iteration parallelism in
/// `importance.rs`) run on this pool.
pub fn run_simulation(
    pedigree: Pedigree,
    marker_set: &crate::MarkerSet,
    params: SimulationParameters,
    progress_tx: Option<std::sync::mpsc::Sender<ProgressEvent>>,
    cancel_flag: Option<Arc<AtomicBool>>,
) -> Result<SimulationResult> {
    let n_threads = params.number_of_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|e| crate::MatchyError::SimulationFailed(e.to_string()))?;
    pool.install(|| run_simulation_impl(pedigree, marker_set, &params, progress_tx, cancel_flag))
}

fn run_simulation_impl(
    mut pedigree: Pedigree,
    marker_set: &crate::MarkerSet,
    params: &SimulationParameters,
    progress_tx: Option<std::sync::mpsc::Sender<ProgressEvent>>,
    cancel_flag: Option<Arc<AtomicBool>>,
) -> Result<SimulationResult> {
    validate_dag(&pedigree)?;

    // -----------------------------------------------------------------------
    // Phase 0: Identify/install comparison haplotype
    // -----------------------------------------------------------------------

    let suspect_haplotype: Option<crate::Haplotype> = if params.trace_mode {
        // Trace mode: add the trace haplotype as a child of the first known
        // individual, exactly as Python's add_trace does.
        let trace = pedigree.trace_haplotype.clone();
        match trace {
            Some(t) => {
                let trace_name = pedigree.add_trace(t.clone());
                if trace_name.is_none() {
                    tracing::warn!("Trace mode: no known individual to attach trace to");
                }
                Some(t)
            }
            None => {
                tracing::warn!("Trace mode enabled but no TRACE profile found in haplotypes JSON");
                None
            }
        }
    } else if let Some(ref suspect_name) = params.suspect {
        let ind = pedigree
            .get_individual_by_name(suspect_name)
            .ok_or_else(|| crate::MatchyError::IndividualNotFound(suspect_name.clone()))?;
        if ind.haplotype.alleles.is_empty() {
            return Err(crate::MatchyError::SimulationFailed(format!(
                "Suspect '{}' has no typed haplotype in the JSON file",
                suspect_name
            )));
        }
        let hap = ind.haplotype.clone();
        let suspect_id = ind.id.clone();
        // Mark suspect as Suspect class — required for remove_irrelevant_individuals
        // and calculate_picking_probabilities to identify the anchor individual.
        if let Some(ind_mut) = pedigree.get_individual_by_id_mut(&suspect_id) {
            ind_mut.haplotype_class = HaplotypeClass::Suspect;
        }
        Some(hap)
    } else {
        None
    };

    // Mark excluded individuals based on params.exclude (list of names).
    for name in &params.exclude {
        if let Some(ind) = pedigree.get_individual_by_name_mut(name) {
            ind.exclude = true;
        }
    }

    // -----------------------------------------------------------------------
    // Phase 1: Prepare pedigrees (extend, prune, reroot)
    // Mirrors Python simulation.py:697–773.
    // -----------------------------------------------------------------------

    // Clone before any structural changes for the extended (outside) pedigree.
    let mut ext_ped = pedigree.clone();

    // Extend the clone — must be called on the original (un-rerooted) DAG so
    // that BFS depth is measured from the true root.  Mirrors Python:700.
    let last_child_name = ext_ped.extend_pedigree();

    // Prune irrelevant known individuals and get the root name to reroot at.
    // For the inside pedigree.
    let root_name = pedigree.remove_irrelevant_individuals(true, None);
    let ext_root_name = ext_ped.remove_irrelevant_individuals(false, Some(&last_child_name));

    // Reroot both pedigrees.
    if let Some(ref name) = root_name {
        if let Some(ind) = pedigree.get_individual_by_name(name) {
            let rid = ind.id.clone();
            pedigree.reroot(&rid);
        }
    }
    if let Some(ref name) = ext_root_name {
        if let Some(ind) = ext_ped.get_individual_by_name(name) {
            let rid = ind.id.clone();
            ext_ped.reroot(&rid);
        }
    }

    // -----------------------------------------------------------------------
    // Phase 2: Picking probabilities
    // -----------------------------------------------------------------------

    if suspect_haplotype.is_some() {
        pedigree.calculate_picking_probabilities(params.trace_mode);

        // Extended pedigree: only the last added child gets probability 1.
        ext_ped.picking_probabilities.clear();
        if let Some(last_child) = ext_ped.get_individual_by_name(&last_child_name) {
            let lcid = last_child.id.clone();
            ext_ped.picking_probabilities.insert(lcid, Decimal::ONE);
        }
    }

    // -----------------------------------------------------------------------
    // Phase 3: BFS setup for original pedigree
    // -----------------------------------------------------------------------

    let roots = pedigree.roots();
    let root_id = roots
        .first()
        .ok_or_else(|| crate::MatchyError::InvalidPedigree("No root individual".into()))?
        .to_string();

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

    // -----------------------------------------------------------------------
    // Step 1: Average pedigree probability
    // -----------------------------------------------------------------------

    let prog_ref = progress_tx.as_ref();
    let cancel = cancel_flag.as_deref();
    let sim_start = Instant::now();

    let t_ped = Instant::now();
    let pedigree_prob_trial = if known_count <= 1 {
        tracing::info!("Only one known individual — pedigree probability is 1");
        let mut trial = convergence::EnsembleTrial::new(1);
        trial.converged = true;
        trial.grand_mean = Some(Decimal::ONE);
        trial
    } else {
        run_ensemble_pedigree_probability(
            &pedigree,
            &root_id,
            marker_set,
            params,
            prog_ref,
            adaptive.as_mut(),
            cancel,
            SimulationStage::PedigreeProbability,
        )?
    };
    let pedigree_elapsed = t_ped.elapsed();
    let pedigree_stats = make_stage_stats(&pedigree_prob_trial, params.batch_length, pedigree_elapsed);

    let avg_pedigree_probability = pedigree_prob_trial.grand_mean.unwrap_or(Decimal::ZERO);

    // -----------------------------------------------------------------------
    // Step 2: Inside-pedigree match probabilities
    // -----------------------------------------------------------------------

    let mut inside_stats: Option<StageStats> = None;
    let (inside_match_probabilities, per_individual_probabilities) =
        if let Some(ref suspect_hap) = suspect_haplotype {
            if params.skip_inside {
                (None, None)
            } else {
                let t_inside = Instant::now();
                let trial = run_ensemble_matching_haplotypes(
                    &pedigree,
                    &root_id,
                    suspect_hap,
                    avg_pedigree_probability,
                    marker_set,
                    params,
                    false, // inside
                    prog_ref,
                    adaptive.as_mut(),
                    cancel,
                )?;
                inside_stats = Some(make_stage_stats(&trial, params.batch_length, t_inside.elapsed()));
                let probs = aggregate_match_counts(&trial.model_results);
                let per_indiv_raw = aggregate_per_individual(&trial.model_results);
                let per_indiv_named: std::collections::HashMap<String, Decimal> = per_indiv_raw
                    .into_iter()
                    .filter_map(|(id, prob)| {
                        pedigree.get_individual_by_id(&id).map(|ind| (ind.name.clone(), prob))
                    })
                    .collect();
                (
                    Some(MatchProbabilities {
                        probabilities: probs,
                        average_pedigree_probability: avg_pedigree_probability,
                    }),
                    if per_indiv_named.is_empty() { None } else { Some(per_indiv_named) },
                )
            }
        } else {
            (None, None)
        };

    // -----------------------------------------------------------------------
    // Step 3: Outside-pedigree match probability
    //
    // Uses the separately prepared extended pedigree.
    // Mirrors Python simulation.py:873–927.
    // -----------------------------------------------------------------------

    let mut extended_pedigree_stats: Option<StageStats> = None;
    let mut outside_stats: Option<StageStats> = None;
    let outside_match_probability = if let Some(ref suspect_hap) = suspect_haplotype {
        if params.skip_outside {
            None
        } else {
            let ext_roots = ext_ped.roots();
            let ext_root_id = ext_roots
                .first()
                .ok_or_else(|| crate::MatchyError::InvalidPedigree("Extended pedigree has no root".into()))?
                .to_string();

            let ext_known_count = {
                let ext_layers = bfs_layers(&ext_ped, &ext_root_id).unwrap_or_default();
                ext_layers.into_iter().flatten().filter(|id| {
                    ext_ped.get_individual_by_id(id)
                        .map(|i| i.haplotype_class != HaplotypeClass::Unknown)
                        .unwrap_or(false)
                }).count()
            };
            let t_ext_ped = Instant::now();
            let ext_avg_ped_prob = if ext_known_count <= 1 {
                Decimal::ONE
            } else {
                let ext_trial = run_ensemble_pedigree_probability(
                    &ext_ped, &ext_root_id, marker_set, params, prog_ref, adaptive.as_mut(), cancel,
                    SimulationStage::ExtendedPedigreeProbability,
                )?;
                extended_pedigree_stats = Some(make_stage_stats(&ext_trial, params.batch_length, t_ext_ped.elapsed()));
                ext_trial.grand_mean.unwrap_or(Decimal::ZERO)
            };

            let t_outside = Instant::now();
            let trial = run_ensemble_matching_haplotypes(
                &ext_ped,
                &ext_root_id,
                suspect_hap,
                ext_avg_ped_prob,
                marker_set,
                params,
                true, // is_outside
                prog_ref,
                adaptive.as_mut(),
                cancel,
            )?;
            outside_stats = Some(make_stage_stats(&trial, params.batch_length, t_outside.elapsed()));
            trial.grand_mean
        }
    } else {
        None
    };

    Ok(SimulationResult {
        parameters: params.clone(),
        inside_match_probabilities,
        outside_match_probability,
        per_individual_probabilities,
        trials: pedigree_prob_trial.trial_nr,
        converged: pedigree_prob_trial.converged,
        pedigree_stats,
        inside_stats,
        extended_pedigree_stats,
        outside_stats,
        total_runtime_secs: sim_start.elapsed().as_secs_f64(),
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
    ExtendedPedigreeProbability,
    InsideMatchProbability,
    OutsideMatchProbability,
}
