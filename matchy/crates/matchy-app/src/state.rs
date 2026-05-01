use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicBool, Ordering};

/// Global application state managed by Tauri.
#[derive(Default)]
pub struct AppState {
    /// Current simulation handle, if one is running.
    pub simulation: Mutex<Option<SimulationHandle>>,
}

/// Handle to a running simulation.
pub struct SimulationHandle {
    pub cancel_flag: Arc<AtomicBool>,
}

impl SimulationHandle {
    pub fn new() -> (Self, Arc<AtomicBool>) {
        let flag = Arc::new(AtomicBool::new(false));
        (Self { cancel_flag: flag.clone() }, flag)
    }

    pub fn cancel(&self) {
        self.cancel_flag.store(true, Ordering::Relaxed);
    }
}
