pub mod html;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum ReportError {
    #[error("Template error: {0}")]
    Template(#[from] minijinja::Error),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

pub type Result<T> = std::result::Result<T, ReportError>;
