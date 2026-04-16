pub mod models;
pub mod graph;
pub mod hungarian;
pub mod simulation;

pub use models::*;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum MatchyError {
    #[error("Invalid pedigree: {0}")]
    InvalidPedigree(String),

    #[error("Individual not found: {0}")]
    IndividualNotFound(String),

    #[error("Marker not found: {0}")]
    MarkerNotFound(String),

    #[error("Simulation failed: {0}")]
    SimulationFailed(String),

    #[error("Convergence failed: {0}")]
    ConvergenceFailed(String),

    #[error("Decimal arithmetic error: {0}")]
    DecimalError(String),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

pub type Result<T> = std::result::Result<T, MatchyError>;
