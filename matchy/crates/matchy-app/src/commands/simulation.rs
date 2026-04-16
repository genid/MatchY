use matchy_core::SimulationParameters;
use serde::{Deserialize, Serialize};
use tauri::{AppHandle, State};

use crate::state::AppState;

/// Parameters passed from the frontend to start a simulation.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulationRequest {
    pub pedigree_tgf: String,
    pub haplotypes_json: String,
    pub marker_set_name: Option<String>,
    pub marker_set_csv: Option<String>,
    pub suspect: Option<String>,
    pub exclude: Vec<String>,
    pub two_step_mutation_fraction: f64,
    pub batch_length: u64,
    pub convergence_criterion: f64,
    pub bias: Option<f64>,
    pub number_of_threads: usize,
    pub skip_inside: bool,
    pub skip_outside: bool,
    pub trace_mode: bool,
    pub adaptive_bias: bool,
    pub simulation_name: String,
    pub user_name: String,
}

/// Result returned to the frontend when simulation completes.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulationResponse {
    pub success: bool,
    pub error: Option<String>,
    // TODO: add structured result fields (match probabilities, etc.) in Phase 4
}

/// Start a simulation. Progress events are emitted via Tauri events.
#[tauri::command]
pub async fn run_simulation(
    app_handle: AppHandle,
    state: State<'_, AppState>,
    request: SimulationRequest,
) -> Result<SimulationResponse, String> {
    use crate::state::SimulationHandle;

    // Cancel any running simulation
    {
        let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
        if let Some(handle) = sim_lock.take() {
            handle.cancel();
        }
    }

    let (handle, _cancel_rx) = SimulationHandle::new();
    {
        let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
        *sim_lock = Some(handle);
    }

    // TODO: Phase 4 — wire up the actual simulation
    // 1. Parse pedigree from request.pedigree_tgf
    // 2. Load marker set
    // 3. Parse haplotypes
    // 4. Build SimulationParameters
    // 5. tokio::task::spawn_blocking(move || matchy_core::simulation::run_simulation(...))
    // 6. Stream ProgressEvents via app_handle.emit("simulation-progress", event)

    let _ = (app_handle, request);

    Ok(SimulationResponse {
        success: false,
        error: Some("Simulation engine not yet implemented (Phase 2/4)".into()),
    })
}

/// Cancel the currently running simulation.
#[tauri::command]
pub async fn cancel_simulation(state: State<'_, AppState>) -> Result<(), String> {
    let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
    if let Some(handle) = sim_lock.take() {
        handle.cancel();
    }
    Ok(())
}
