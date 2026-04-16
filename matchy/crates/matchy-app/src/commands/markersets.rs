use serde::Serialize;

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct MarkerInfo {
    pub name: String,
    pub mutation_rate: f64,
    pub number_of_copies: Option<u32>,
}

/// List all available built-in kit names.
#[tauri::command]
pub fn list_kits() -> Result<Vec<String>, String> {
    matchy_io::kits::kit_names().map_err(|e| e.to_string())
}

/// Load a built-in kit by name and return its markers.
#[tauri::command]
pub fn load_kit(name: String) -> Result<Vec<MarkerInfo>, String> {
    let marker_set = matchy_io::kits::load_kit(&name)
        .map_err(|e| e.to_string())?
        .ok_or_else(|| format!("Kit '{}' not found", name))?;

    Ok(marker_set
        .markers
        .iter()
        .map(|m| MarkerInfo {
            name: m.name.clone(),
            mutation_rate: m.mutation_rate,
            number_of_copies: m.number_of_copies,
        })
        .collect())
}

/// Load a custom marker set from a CSV string.
#[tauri::command]
pub fn load_custom_csv(csv_content: String) -> Result<Vec<MarkerInfo>, String> {
    let marker_set = matchy_io::marker_csv::read_marker_csv_str(&csv_content)
        .map_err(|e| e.to_string())?;

    Ok(marker_set
        .markers
        .iter()
        .map(|m| MarkerInfo {
            name: m.name.clone(),
            mutation_rate: m.mutation_rate,
            number_of_copies: m.number_of_copies,
        })
        .collect())
}
