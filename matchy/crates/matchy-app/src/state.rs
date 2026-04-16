use std::sync::Mutex;
use tokio::sync::broadcast;

use crate::progress::ProgressEvent;

/// Global application state managed by Tauri.
#[derive(Default)]
pub struct AppState {
    /// Current simulation handle, if one is running.
    pub simulation: Mutex<Option<SimulationHandle>>,
}

/// Handle to a running simulation.
pub struct SimulationHandle {
    /// Send on this channel to request cancellation.
    pub cancel_tx: broadcast::Sender<()>,
}

impl SimulationHandle {
    pub fn new() -> (Self, broadcast::Receiver<()>) {
        let (cancel_tx, cancel_rx) = broadcast::channel(1);
        (Self { cancel_tx }, cancel_rx)
    }

    pub fn cancel(&self) {
        let _ = self.cancel_tx.send(());
    }
}
