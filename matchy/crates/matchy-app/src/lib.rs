pub mod commands;
pub mod progress;
pub mod state;

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .manage(state::AppState::default())
        .invoke_handler(tauri::generate_handler![
            commands::simulation::run_simulation,
            commands::simulation::cancel_simulation,
            commands::simulation::get_cpu_count,
            commands::pedigree::parse_tgf,
            commands::pedigree::parse_ped,
            commands::pedigree::export_tgf,
            commands::pedigree::validate_pedigree,
            commands::pedigree::build_extended_pedigree,
            commands::haplotypes::parse_haplotypes_json,
            commands::haplotypes::export_haplotypes_json,
            commands::markersets::list_kits,
            commands::markersets::list_all_markers,
            commands::markersets::load_kit,
            commands::markersets::load_custom_csv,
            commands::report::generate_report,
            commands::report::save_and_open_report,
            commands::report::save_run,
        ])
        .setup(|app| {
            #[cfg(debug_assertions)]
            {
                use tauri::Manager;
                app.get_webview_window("main").unwrap().open_devtools();
            }
            let _ = app;
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running Tauri application");
}
