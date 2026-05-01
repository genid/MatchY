pub mod config;
pub mod haplotypes_json;
pub mod kits;
pub mod marker_csv;
pub mod mutation_rates;
pub mod ped;
pub mod tgf;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum IoError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("JSON parse error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("TOML parse error: {0}")]
    Toml(#[from] toml::de::Error),

    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),

    #[error("Invalid file format: {0}")]
    InvalidFormat(String),

    #[error("Individual '{0}' not found in pedigree")]
    IndividualNotFound(String),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

pub type Result<T> = std::result::Result<T, IoError>;
