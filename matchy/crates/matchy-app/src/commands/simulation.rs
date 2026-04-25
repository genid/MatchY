use matchy_core::{SimulationParameters, SimulationResult};
use serde::{Deserialize, Serialize};
use tauri::{AppHandle, Emitter, State};

use crate::progress::ProgressEvent;
use crate::state::AppState;

// ---------------------------------------------------------------------------
// Request / Response types
// ---------------------------------------------------------------------------

/// Parameters sent from the frontend to start a simulation.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulationRequest {
    /// TGF content as a string (loaded/edited in the Pedigree Builder).
    pub pedigree_tgf: String,
    /// Haplotypes JSON content as a string.
    pub haplotypes_json: String,
    /// Embedded kit name (e.g. "RMplex"), or None if marker_set_csv is set.
    pub marker_set_name: Option<String>,
    /// Custom marker set CSV content, or None if marker_set_name is set.
    pub marker_set_csv: Option<String>,
    pub suspect: Option<String>,
    #[serde(default)]
    pub exclude: Vec<String>,
    #[serde(default = "default_two_step")]
    pub two_step_mutation_fraction: f64,
    #[serde(default = "default_batch")]
    pub batch_length: u64,
    #[serde(default = "default_convergence")]
    pub convergence_criterion: f64,
    pub bias: Option<f64>,
    #[serde(default = "default_threads")]
    pub number_of_threads: usize,
    #[serde(default)]
    pub skip_inside: bool,
    #[serde(default)]
    pub skip_outside: bool,
    #[serde(default)]
    pub trace_mode: bool,
    #[serde(default)]
    pub adaptive_bias: bool,
    #[serde(default = "default_sim_name")]
    pub simulation_name: String,
    #[serde(default)]
    pub user_name: String,
    pub seed: Option<u64>,
}

fn default_two_step() -> f64 { 0.03 }
fn default_batch() -> u64 { 10_000 }
fn default_convergence() -> f64 { 0.02 }
fn default_threads() -> usize { 3 }
fn default_sim_name() -> String { "simulation".into() }

/// Returned to the frontend when simulation completes or fails.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SimulationResponse {
    pub success: bool,
    pub error: Option<String>,
    /// Structured result — present on success.
    pub result: Option<SimulationResult>,
}

// ---------------------------------------------------------------------------
// run_simulation
// ---------------------------------------------------------------------------

/// Start a simulation.
///
/// Progress events are emitted via Tauri event `"simulation-progress"` as the
/// 3-model ensemble runs. The frontend can `listen("simulation-progress", ...)`
/// to stream live convergence data.
///
/// Returns the full `SimulationResult` on success, or an error string on failure.
#[tauri::command]
pub async fn run_simulation(
    app_handle: AppHandle,
    state: State<'_, AppState>,
    request: SimulationRequest,
) -> Result<SimulationResponse, String> {
    use crate::state::SimulationHandle;

    // Cancel any previously running simulation.
    {
        let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
        if let Some(handle) = sim_lock.take() {
            handle.cancel();
        }
    }

    let (handle, cancel_flag) = SimulationHandle::new();
    {
        let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
        *sim_lock = Some(handle);
    }

    // -----------------------------------------------------------------------
    // Parse inputs on the current async task (fast, no I/O blocking)
    // -----------------------------------------------------------------------

    // 1. Pedigree
    let mut pedigree = matchy_io::tgf::read_tgf_str(&request.pedigree_tgf)
        .map_err(|e| format!("Failed to parse pedigree: {e}"))?;

    // 2. Marker set
    let mut marker_set = if let Some(ref kit_name) = request.marker_set_name {
        matchy_io::kits::load_kit(kit_name)
            .map_err(|e| format!("Failed to load kit: {e}"))?
            .ok_or_else(|| format!("Kit '{}' not found", kit_name))?
    } else if let Some(ref csv) = request.marker_set_csv {
        matchy_io::marker_csv::read_marker_csv(std::io::Cursor::new(csv.as_bytes()))
            .map_err(|e| format!("Failed to parse marker set CSV: {e}"))?
    } else {
        return Ok(SimulationResponse {
            success: false,
            error: Some("Either markerSetName or markerSetCsv must be provided".into()),
            result: None,
        });
    };

    // 3. Haplotypes
    matchy_io::haplotypes_json::read_haplotypes(
        std::io::Cursor::new(request.haplotypes_json.as_bytes()),
        &mut pedigree,
        &mut marker_set,
    )
    .map_err(|e| format!("Failed to parse haplotypes: {e}"))?;

    // 4. SimulationParameters
    let params = SimulationParameters {
        two_step_mutation_fraction: request.two_step_mutation_fraction,
        batch_length: request.batch_length,
        convergence_criterion: request.convergence_criterion,
        bias: request.bias,
        number_of_threads: request.number_of_threads,
        suspect: request.suspect,
        exclude: request.exclude,
        skip_inside: request.skip_inside,
        skip_outside: request.skip_outside,
        trace_mode: request.trace_mode,
        adaptive_bias: request.adaptive_bias,
        simulation_name: request.simulation_name,
        user_name: request.user_name,
        results_path: std::path::PathBuf::from("."),
        seed: request.seed,
    };

    // -----------------------------------------------------------------------
    // Run simulation on a dedicated blocking thread.
    // Progress events are streamed via mpsc → Tauri event.
    // -----------------------------------------------------------------------

    let (progress_tx, progress_rx) = std::sync::mpsc::channel::<matchy_core::simulation::ProgressEvent>();

    // Spawn a task that drains the progress channel and emits Tauri events.
    let app_handle_emit = app_handle.clone();
    let emitter_handle = tokio::task::spawn_blocking(move || {
        while let Ok(core_event) = progress_rx.recv() {
            let fe_event: ProgressEvent = core_event.into();
            let _ = app_handle_emit.emit("simulation-progress", &fe_event);
        }
    });

    // Run the simulation on a separate blocking thread (CPU-bound work must
    // not run on the tokio executor).
    let sim_result = tokio::task::spawn_blocking(move || {
        matchy_core::simulation::run_simulation(pedigree, &marker_set, params, Some(progress_tx), Some(cancel_flag))
    })
    .await
    .map_err(|e| format!("Simulation task panicked: {e}"))?;

    // Wait for the emitter to flush all remaining events.
    let _ = emitter_handle.await;

    match sim_result {
        Ok(result) => {
            // Emit a final "simulation-complete" event so the frontend can transition state.
            let _ = app_handle.emit("simulation-complete", &result);
            Ok(SimulationResponse {
                success: true,
                error: None,
                result: Some(result),
            })
        }
        Err(matchy_core::MatchyError::Cancelled) => Ok(SimulationResponse {
            success: false,
            error: Some("cancelled".into()),
            result: None,
        }),
        Err(e) => Ok(SimulationResponse {
            success: false,
            error: Some(e.to_string()),
            result: None,
        }),
    }
}

// ---------------------------------------------------------------------------
// cancel_simulation
// ---------------------------------------------------------------------------

/// Return the number of logical CPU cores available on this machine.
#[tauri::command]
pub fn get_cpu_count() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

/// Cancel the currently running simulation.
///
/// Note: cancellation is cooperative — the simulation will stop cleanly after
/// completing the current batch, not mid-computation.
#[tauri::command]
pub async fn cancel_simulation(state: State<'_, AppState>) -> Result<(), String> {
    let mut sim_lock = state.simulation.lock().map_err(|e| e.to_string())?;
    if let Some(handle) = sim_lock.take() {
        handle.cancel();
    }
    Ok(())
}
