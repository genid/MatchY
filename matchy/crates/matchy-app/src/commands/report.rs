use std::collections::HashMap;

#[tauri::command]
pub fn generate_report(
    result_json: String,
    pedigree_image: Option<String>,
    chart_images: HashMap<String, String>,
    pedigree_json: Option<String>,
    haplotypes_json: Option<String>,
    markers_json: Option<String>,
    report_date: Option<String>,
    progress_events_json: Option<String>,
) -> Result<String, String> {
    let result: matchy_core::SimulationResult =
        serde_json::from_str(&result_json).map_err(|e| e.to_string())?;

    let html = if result.parameters.trace_mode {
        matchy_report::html::render_trace_report(
            &result,
            pedigree_image.as_deref(),
            &chart_images,
            pedigree_json.as_deref(),
            haplotypes_json.as_deref(),
            markers_json.as_deref(),
            report_date.as_deref(),
            progress_events_json.as_deref(),
        )
    } else {
        matchy_report::html::render_report(
            &result,
            pedigree_image.as_deref(),
            &chart_images,
            pedigree_json.as_deref(),
            haplotypes_json.as_deref(),
            markers_json.as_deref(),
            report_date.as_deref(),
            progress_events_json.as_deref(),
        )
    };

    html.map_err(|e| e.to_string())
}

#[tauri::command]
pub fn save_and_open_report(
    html: String,
    simulation_name: String,
) -> Result<String, String> {
    use std::io::Write;

    let temp_dir = std::env::temp_dir();
    let safe_name: String = simulation_name
        .chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect();
    let file_name = format!("matchy_report_{}.html", safe_name);
    let file_path = temp_dir.join(&file_name);

    {
        let mut file = std::fs::File::create(&file_path)
            .map_err(|e| format!("Failed to create report file: {}", e))?;
        file.write_all(html.as_bytes())
            .map_err(|e| format!("Failed to write report: {}", e))?;
    }

    let path_str = file_path.to_string_lossy().to_string();

    #[cfg(target_os = "windows")]
    {
        std::process::Command::new("cmd")
            .args(["/c", "start", "", &path_str])
            .spawn()
            .map_err(|e| format!("Failed to open report: {}", e))?;
    }
    #[cfg(target_os = "macos")]
    {
        std::process::Command::new("open")
            .arg(&path_str)
            .spawn()
            .map_err(|e| format!("Failed to open report: {}", e))?;
    }
    #[cfg(target_os = "linux")]
    {
        std::process::Command::new("xdg-open")
            .arg(&path_str)
            .spawn()
            .map_err(|e| format!("Failed to open report: {}", e))?;
    }

    Ok(path_str)
}

fn unix_secs_to_datetime_str(secs: u64) -> String {
    let sec = secs % 60;
    let min = (secs / 60) % 60;
    let hour = (secs / 3600) % 24;
    let days = (secs / 86400) as i64;
    // Civil calendar algorithm (Howard Hinnant)
    let z = days + 719468;
    let era = if z >= 0 { z / 146097 } else { (z - 146096) / 146097 };
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let mon = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if mon <= 2 { y + 1 } else { y };
    format!("{:04}-{:02}-{:02}_{:02}-{:02}-{:02}", year, mon, d, hour, min, sec)
}

#[tauri::command]
pub fn save_run(
    folder: String,
    html: String,
    result_json: String,
    simulation_name: String,
) -> Result<String, String> {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let datetime = unix_secs_to_datetime_str(secs);
    let safe_name: String = simulation_name
        .chars()
        .map(|c| if c.is_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect();
    let dir_name = format!("{}_{}", datetime, safe_name);

    let run_dir = std::path::Path::new(&folder).join(&dir_name);
    std::fs::create_dir_all(&run_dir)
        .map_err(|e| format!("Failed to create run folder: {}", e))?;

    std::fs::write(run_dir.join("report.html"), html.as_bytes())
        .map_err(|e| format!("Failed to write report.html: {}", e))?;

    std::fs::write(run_dir.join("result.json"), result_json.as_bytes())
        .map_err(|e| format!("Failed to write result.json: {}", e))?;

    Ok(run_dir.to_string_lossy().to_string())
}
