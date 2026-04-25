/// TOML configuration loader with legacy INI support.
///
/// New format: config.toml (see plan for schema).
/// Legacy: config.ini — detected automatically and converted.
use crate::{IoError, Result};
use matchy_core::SimulationParameters;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

// ---------------------------------------------------------------------------
// TOML schema
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize, Deserialize)]
pub struct Config {
    pub simulation: SimulationSection,
    pub files: FilesSection,
    pub parameters: ParametersSection,
    #[serde(rename = "mode")]
    pub mode: ModeSection,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct SimulationSection {
    pub name: String,
    #[serde(default)]
    pub user: String,
    #[serde(default = "default_results_path")]
    pub results_path: PathBuf,
}

fn default_results_path() -> PathBuf {
    PathBuf::from("./results")
}

#[derive(Debug, Serialize, Deserialize)]
pub struct FilesSection {
    pub pedigree: PathBuf,
    pub known_haplotypes: PathBuf,
    /// Kit name (e.g. "Yfiler plus") OR omit and use marker_set_file
    #[serde(default)]
    pub marker_set: Option<String>,
    /// Path to a custom CSV marker set (used when marker_set is None)
    #[serde(default)]
    pub marker_set_file: Option<PathBuf>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ParametersSection {
    #[serde(default = "default_two_step")]
    pub two_step_mutation_fraction: f64,
    #[serde(default = "default_batch")]
    pub batch_length: u64,
    #[serde(default = "default_convergence")]
    pub convergence_criterion: f64,
    #[serde(default)]
    pub number_of_threads: Option<usize>,
    #[serde(default)]
    pub bias: Option<f64>,
    #[serde(default)]
    pub seed: Option<u64>,
}

fn default_two_step() -> f64 { 0.03 }
fn default_batch() -> u64 { 10_000 }
fn default_convergence() -> f64 { 0.02 }

#[derive(Debug, Serialize, Deserialize)]
pub struct ModeSection {
    #[serde(default)]
    pub suspect: Option<String>,
    #[serde(default)]
    pub exclude: Vec<String>,
    #[serde(default)]
    pub skip_inside: bool,
    #[serde(default)]
    pub skip_outside: bool,
    #[serde(default)]
    pub trace_mode: bool,
    #[serde(default)]
    pub adaptive_bias: bool,
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

/// Load a config file. Detects TOML vs legacy INI by extension.
pub fn load_config(path: &Path) -> Result<Config> {
    let content = std::fs::read_to_string(path)?;
    match path.extension().and_then(|e| e.to_str()) {
        Some("toml") => load_toml_str(&content),
        Some("ini") | Some("cfg") => load_ini_legacy(&content),
        _ => {
            // Try TOML first, then INI
            load_toml_str(&content).or_else(|_| load_ini_legacy(&content))
        }
    }
}

/// Parse a TOML config string.
pub fn load_toml_str(s: &str) -> Result<Config> {
    toml::from_str(s).map_err(IoError::Toml)
}

/// Convert a legacy INI config to the Config struct.
/// Supports the original [pedigree] section layout.
pub fn load_ini_legacy(content: &str) -> Result<Config> {
    use std::collections::HashMap;

    // Simple key=value parser (no section headers — original config has [pedigree])
    let mut values: HashMap<String, String> = HashMap::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with('[') {
            continue;
        }
        if let Some(eq_pos) = line.find('=') {
            let key = line[..eq_pos].trim().to_string();
            let value = line[eq_pos + 1..].trim().to_string();
            values.insert(key, value);
        }
    }

    let get = |key: &str| -> Option<String> { values.get(key).cloned() };

    let config = Config {
        simulation: SimulationSection {
            name: get("simulation_name").unwrap_or_else(|| "simulation".into()),
            user: get("user_name").unwrap_or_default(),
            results_path: get("results_path")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("./results")),
        },
        files: FilesSection {
            pedigree: get("path")
                .ok_or_else(|| IoError::InvalidFormat("Missing 'path' in INI config".into()))?
                .into(),
            known_haplotypes: get("known_haplotypes")
                .ok_or_else(|| {
                    IoError::InvalidFormat("Missing 'known_haplotypes' in INI config".into())
                })?
                .into(),
            marker_set: get("marker_set"),
            marker_set_file: None,
        },
        parameters: ParametersSection {
            two_step_mutation_fraction: get("two_step_mutation_fraction")
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.03),
            batch_length: get("batch_length")
                .and_then(|v| v.parse().ok())
                .unwrap_or(10_000),
            convergence_criterion: get("convergence_criterion")
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.02),
            number_of_threads: get("number_of_threads").and_then(|v| v.parse().ok()),
            bias: get("bias").and_then(|v| {
                if v == "None" || v.is_empty() {
                    None
                } else {
                    v.parse().ok()
                }
            }),
            seed: get("seed").and_then(|v| v.parse().ok()),
        },
        mode: ModeSection {
            suspect: get("suspect"),
            exclude: get("exclude")
                .map(|v| v.split(',').map(|s| s.trim().to_string()).collect())
                .unwrap_or_default(),
            skip_inside: false,
            skip_outside: false,
            trace_mode: false,
            adaptive_bias: false,
        },
    };

    Ok(config)
}

// ---------------------------------------------------------------------------
// Convert Config → SimulationParameters
// ---------------------------------------------------------------------------

impl From<Config> for SimulationParameters {
    fn from(c: Config) -> Self {
        let threads = c
            .parameters
            .number_of_threads
            .unwrap_or_else(|| num_cpus());
        SimulationParameters {
            two_step_mutation_fraction: c.parameters.two_step_mutation_fraction,
            batch_length: c.parameters.batch_length,
            convergence_criterion: c.parameters.convergence_criterion,
            bias: c.parameters.bias,
            number_of_threads: threads,
            suspect: c.mode.suspect,
            exclude: c.mode.exclude,
            skip_inside: c.mode.skip_inside,
            skip_outside: c.mode.skip_outside,
            trace_mode: c.mode.trace_mode,
            adaptive_bias: c.mode.adaptive_bias,
            simulation_name: c.simulation.name,
            user_name: c.simulation.user,
            results_path: c.simulation.results_path,
            seed: c.parameters.seed,
        }
    }
}

fn num_cpus() -> usize {
    // Use rayon's default thread count (== logical CPU count)
    rayon::current_num_threads()
}
