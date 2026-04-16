// Prevents additional console window on Windows in release mode.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

fn main() {
    // If command-line arguments are present (besides the binary name and any
    // Tauri-internal flags), delegate to the CLI handler.
    let args: Vec<String> = std::env::args().collect();
    let is_cli = args.len() > 1
        && !args[1].starts_with("--tauri")
        && !args[1].starts_with("tauri");

    if is_cli {
        // Re-invoke as the CLI binary by delegating to matchy-cli's logic.
        // For now, print a placeholder message.
        // TODO: integrate matchy-cli::run(args) here once Phase 2 is complete.
        eprintln!("CLI mode detected. Use the matchy binary for headless operation.");
        std::process::exit(1);
    }

    matchy_app_lib::run();
}
