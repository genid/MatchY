pub mod commands;
pub mod progress;
pub mod state;

use tauri::Manager;

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .manage(state::AppState::default())
        .invoke_handler(tauri::generate_handler![
            commands::simulation::run_simulation,
            commands::simulation::cancel_simulation,
            commands::pedigree::parse_tgf,
            commands::pedigree::parse_ped,
            commands::pedigree::export_tgf,
            commands::pedigree::validate_pedigree,
            commands::haplotypes::parse_haplotypes_json,
            commands::haplotypes::export_haplotypes_json,
            commands::markersets::list_kits,
            commands::markersets::load_kit,
            commands::markersets::load_custom_csv,
            commands::report::generate_report,
        ])
        .setup(|app| {
            #[cfg(debug_assertions)]
            {
                let window = app.get_webview_window("main").unwrap();
                window.open_devtools();
            }
            Ok(())
        })
        .run(tauri::generate_context!())
        .expect("error while running Tauri application");
}
