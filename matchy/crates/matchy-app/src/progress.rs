use serde::Serialize;
use matchy_core::simulation::{ProgressEvent as CoreProgressEvent, SimulationStage};

/// Progress event emitted to the frontend via Tauri's event system.
/// Mirrors matchy_core::simulation::ProgressEvent with JSON-serialisable types.
#[derive(Debug, Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProgressEvent {
    pub trial: u32,
    pub model: u8,
    pub iteration: u64,
    pub current_mean: String,
    pub stage: String,
    pub converged: bool,
}

impl From<CoreProgressEvent> for ProgressEvent {
    fn from(e: CoreProgressEvent) -> Self {
        let stage = match e.stage {
            SimulationStage::PedigreeProbability => "pedigree_probability",
            SimulationStage::ExtendedPedigreeProbability => "extended_pedigree_probability",
            SimulationStage::InsideMatchProbability => "inside_match_probability",
            SimulationStage::OutsideMatchProbability => "outside_match_probability",
        };
        Self {
            trial: e.trial,
            model: e.model,
            iteration: e.iteration,
            current_mean: e.current_mean,
            stage: stage.to_string(),
            converged: e.converged,
        }
    }
}
