/// Embedded Y-STR kit definitions.
///
/// kits.json format: { "KitName": ["MarkerName1", "MarkerName2", ...] }
/// Mutation rates are looked up from mutation_rates.csv.
use crate::mutation_rates::load_mutation_rates;
use crate::Result;
use matchy_core::{Marker, MarkerSet};
use std::collections::HashMap;

// Embed at compile time for a self-contained binary.
const KITS_JSON: &str = include_str!("../assets/kits.json");

/// Load all bundled kits. Returns a map: kit_name → MarkerSet.
pub fn load_all_kits() -> Result<HashMap<String, MarkerSet>> {
    let kits: HashMap<String, Vec<String>> = serde_json::from_str(KITS_JSON)?;
    let rates = load_mutation_rates()?;
    let mut result = HashMap::new();

    for (kit_name, marker_names) in kits {
        let mut ms = MarkerSet::new();
        for marker_name in &marker_names {
            let mutation_rate = rates.get(marker_name.as_str()).copied().unwrap_or_else(|| {
                tracing::warn!(
                    "Marker '{}' not found in mutation_rates.csv — using rate 0.001",
                    marker_name
                );
                0.001
            });
            ms.add_marker(Marker::new(marker_name.as_str(), mutation_rate));
        }
        result.insert(kit_name, ms);
    }
    Ok(result)
}

/// Load a single kit by name.
pub fn load_kit(name: &str) -> Result<Option<MarkerSet>> {
    let mut kits = load_all_kits()?;
    Ok(kits.remove(name))
}

/// List available kit names (sorted).
pub fn kit_names() -> Result<Vec<String>> {
    let kits: HashMap<String, Vec<String>> = serde_json::from_str(KITS_JSON)?;
    let mut names: Vec<String> = kits.into_keys().collect();
    names.sort();
    Ok(names)
}
