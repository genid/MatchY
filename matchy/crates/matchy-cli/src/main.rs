use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::{bail, Context, Result};
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
    if args.skip_inside   { config.mode.skip_inside   = true; }
    if args.skip_outside  { config.mode.skip_outside  = true; }
    if args.trace_mode    { config.mode.trace_mode    = true; }
    if args.adaptive_bias { config.mode.adaptive_bias = true; }

    // Resolve file paths relative to the config file's directory
    let config_dir = args.config_path
        .parent()
        .unwrap_or_else(|| std::path::Path::new("."))
        .to_path_buf();

    let pedigree_path = resolve_path(&config_dir, &config.files.pedigree);
    let haplotypes_path = resolve_path(&config_dir, &config.files.known_haplotypes);
    let results_path = resolve_path(&config_dir, &config.simulation.results_path);

    // Load marker set
    info!("Loading marker set");
    let mut marker_set = if let Some(ref kit_name) = config.files.marker_set {
        matchy_io::kits::load_kit(kit_name)
            .context("Failed to load embedded kit")?
            .ok_or_else(|| anyhow::anyhow!("Kit '{}' not found. Run 'matchy --list-kits' to see available kits.", kit_name))?
    } else if let Some(ref csv_path) = config.files.marker_set_file {
        let csv_path = resolve_path(&config_dir, csv_path);
        let f = std::fs::File::open(&csv_path)
            .with_context(|| format!("Cannot open marker set file {:?}", csv_path))?;
        matchy_io::marker_csv::read_marker_csv(std::io::BufReader::new(f))
            .context("Failed to parse marker set CSV")?
    } else {
        bail!("Config must specify either [files] marker_set (kit name) or marker_set_file (CSV path)");
    };

    // Load pedigree
    info!("Loading pedigree from {:?}", pedigree_path);
    let pedigree_ext = pedigree_path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("");

    let mut pedigree = match pedigree_ext {
        "tgf" => {
            let f = std::fs::File::open(&pedigree_path)
                .with_context(|| format!("Cannot open pedigree {:?}", pedigree_path))?;
            matchy_io::tgf::read_tgf(std::io::BufReader::new(f))
                .context("Failed to parse TGF pedigree")?
        }
        "ped" => {
            let f = std::fs::File::open(&pedigree_path)
                .with_context(|| format!("Cannot open pedigree {:?}", pedigree_path))?;
            matchy_io::ped::read_ped(std::io::BufReader::new(f))
                .context("Failed to parse PED pedigree")?
        }
        other => bail!("Unsupported pedigree format '.{}'. Expected .tgf or .ped", other),
    };

    // Load haplotypes and assign to pedigree
    info!("Loading haplotypes from {:?}", haplotypes_path);
    let f = std::fs::File::open(&haplotypes_path)
        .with_context(|| format!("Cannot open haplotypes {:?}", haplotypes_path))?;
    let trace_haplotype =
        matchy_io::haplotypes_json::read_haplotypes(f, &mut pedigree, &mut marker_set)
            .context("Failed to parse haplotypes JSON")?;

    // Build SimulationParameters from config
    let params: matchy_core::SimulationParameters = config.into();

    // If trace mode: use the TRACE haplotype as the suspect
    if params.trace_mode {
        if let Some(ref trace_hap) = trace_haplotype {
            // Add a synthetic "TRACE" individual to the pedigree
            // (treated as an outside individual in the simulation)
            info!("Trace mode: using TRACE haplotype");
            let _ = trace_hap; // Used via pedigree logic in simulation
        } else {
            bail!("Trace mode enabled but no TRACE haplotype found in haplotypes JSON");
        }
    }

    // Create results directory
    std::fs::create_dir_all(&results_path)
        .with_context(|| format!("Cannot create results directory {:?}", results_path))?;

    info!(
        "Simulation '{}' — {} thread(s), batch_length={}, convergence_criterion={}",
        params.simulation_name,
        params.number_of_threads,
        params.batch_length,
        params.convergence_criterion,
    );

    // Collect progress events in a background thread so the report gets convergence charts
    let (progress_tx, progress_rx) =
        std::sync::mpsc::channel::<matchy_core::simulation::ProgressEvent>();
    let progress_collector =
        std::thread::spawn(move || progress_rx.into_iter().collect::<Vec<_>>());

    // Run simulation
    info!("Running simulation '{}'", params.simulation_name);
    let result = matchy_core::simulation::run_simulation(
        pedigree,
        &marker_set,
        params,
        Some(progress_tx),
        None,
    )
    .context("Simulation failed")?;

    let progress_events = progress_collector.join().expect("progress collector panicked");
    let progress_events_json = serde_json::to_string(&progress_events).ok();

    // Print summary to stdout
    print_summary(&result);

    // Render HTML report
    let html = if result.parameters.trace_mode {
        matchy_report::html::render_trace_report(
            &result,
            None,
            &std::collections::HashMap::new(),
            None,
            None,
            None,
            None,
            progress_events_json.as_deref(),
            None,
        )
        .context("Failed to render trace report")?
    } else {
        matchy_report::html::render_report(
            &result,
            None,
            &std::collections::HashMap::new(),
            None,
            None,
            None,
            None,
            progress_events_json.as_deref(),
            None,
        )
        .context("Failed to render report")?
    };

    let report_filename = format!("{}_report.html", result.parameters.simulation_name);
    let report_path = results_path.join(&report_filename);
    std::fs::write(&report_path, &html)
        .with_context(|| format!("Cannot write report to {:?}", report_path))?;

    info!("Report written to {:?}", report_path);
    info!("Done.");
    Ok(())
}

/// Resolve a path relative to a base directory (if not already absolute).
fn resolve_path(base: &std::path::Path, path: &std::path::Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        base.join(path)
    }
}

/// Print a human-readable simulation summary.
fn print_summary(result: &matchy_core::SimulationResult) {
    println!("\n=== MatchY Simulation Results ===");
    println!("Simulation : {}", result.parameters.simulation_name);
    println!("Converged  : {} (trials: {})", result.converged, result.trials);

    if let Some(ref inside) = result.inside_match_probabilities {
        println!("\nInside-pedigree match probabilities:");
        println!(
            "  Average pedigree probability: {}",
            inside.average_pedigree_probability
        );
        let p_at_least_one: f64 = inside.probabilities.values()
            .filter_map(|d| f64::try_from(*d).ok())
            .sum();
        println!("  P(at least 1 match) = {:.6}", p_at_least_one);
    }

    if let Some(outside_prob) = result.outside_match_probability {
        println!("\nOutside-pedigree match probability:");
        println!("  P(match outside) = {}", outside_prob);
    }

    if let Some(ref per_ind) = result.per_individual_probabilities {
        println!("\nPer-individual probabilities (trace mode):");
        let mut entries: Vec<_> = per_ind.iter().collect();
        entries.sort_by(|a, b| b.1.partial_cmp(a.1).unwrap_or(std::cmp::Ordering::Equal));
        for (id, prob) in entries {
            println!("  {}: {}", id, prob);
        }
    }

    println!();
}
