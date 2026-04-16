use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;
use tracing::info;

/// MatchY — Y-STR match probability estimator
#[derive(Parser, Debug)]
#[command(name = "matchy", version, about, long_about = None)]
pub struct Args {
    /// Path to configuration file (.toml or legacy .ini)
    #[arg(short = 'c', long, default_value = "config.toml")]
    pub config_path: PathBuf,

    /// Skip inside-pedigree probability calculation
    #[arg(short = 'i', long, default_value_t = false)]
    pub skip_inside: bool,

    /// Skip outside-pedigree probability calculation
    #[arg(short = 'o', long, default_value_t = false)]
    pub skip_outside: bool,

    /// Run in trace mode (identify most likely donor)
    #[arg(short = 't', long, default_value_t = false)]
    pub trace_mode: bool,

    /// Enable adaptive bias adjustment across trials
    #[arg(short = 'a', long, default_value_t = false)]
    pub adaptive_bias: bool,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let args = Args::parse();

    info!("Loading config from {:?}", args.config_path);
    let mut config = matchy_io::config::load_config(&args.config_path)
        .context("Failed to load configuration")?;

    // CLI flags override config file values
    if args.skip_inside  { config.mode.skip_inside  = true; }
    if args.skip_outside { config.mode.skip_outside = true; }
    if args.trace_mode   { config.mode.trace_mode   = true; }
    if args.adaptive_bias { config.mode.adaptive_bias = true; }

    let mut params: matchy_core::SimulationParameters = config.into();

    // Load pedigree
    info!("Loading pedigree from {:?}", params.results_path);
    // TODO: load pedigree from config.files.pedigree

    // Load marker set
    // TODO: load from config.files.marker_set or config.files.marker_set_file

    // Load haplotypes
    // TODO: load from config.files.known_haplotypes

    // Run simulation
    info!("Running simulation '{}'", params.simulation_name);
    // TODO: call matchy_core::simulation::run_simulation(...)

    // Generate report
    // TODO: call matchy_report::html::render_report(...)

    info!("Done. Results written to {:?}", params.results_path);
    Ok(())
}
