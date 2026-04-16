use serde::Deserialize;
use std::collections::HashMap;

/// Generate an HTML report from a completed simulation result.
/// `chart_images` is a map of chart_name → base64-encoded PNG from Chart.js.
/// `pedigree_image` is the base64-encoded PNG from react-flow.
#[tauri::command]
pub fn generate_report(
    result_json: String,
    pedigree_image: Option<String>,
    chart_images: HashMap<String, String>,
) -> Result<String, String> {
    let result: matchy_core::SimulationResult =
        serde_json::from_str(&result_json).map_err(|e| e.to_string())?;

    let html = if result.parameters.trace_mode {
        matchy_report::html::render_trace_report(
            &result,
            pedigree_image.as_deref(),
            &chart_images,
        )
    } else {
        matchy_report::html::render_report(
            &result,
            pedigree_image.as_deref(),
            &chart_images,
        )
    };

    html.map_err(|e| e.to_string())
}
