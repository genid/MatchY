use matchy_core::{MarkerSet, Pedigree};
use serde::{Deserialize, Serialize};
use serde_json::Value;

/// Parse a haplotypes JSON string and return both the updated pedigree data
/// and the trace haplotype (if present).
#[tauri::command]
pub fn parse_haplotypes_json(
    json_content: String,
    pedigree_tgf: String,
    marker_set_name: Option<String>,
    marker_set_csv: Option<String>,
) -> Result<HaplotypesParseResult, String> {
    let mut pedigree =
        matchy_io::tgf::read_tgf_str(&pedigree_tgf).map_err(|e| e.to_string())?;

    let mut marker_set = if let Some(name) = marker_set_name {
        matchy_io::kits::load_kit(&name)
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("Kit '{}' not found", name))?
    } else if let Some(csv) = marker_set_csv {
        matchy_io::marker_csv::read_marker_csv_str(&csv).map_err(|e| e.to_string())?
    } else {
        return Err("Either marker_set_name or marker_set_csv must be provided".into());
    };

    let trace = matchy_io::haplotypes_json::read_haplotypes(
        json_content.as_bytes(),
        &mut pedigree,
        &mut marker_set,
    )
    .map_err(|e| e.to_string())?;

    // Build frontend-friendly haplotype table: individual_name → marker_name → allele_string
    let mut table: std::collections::HashMap<String, std::collections::HashMap<String, String>> =
        std::collections::HashMap::new();

    for individual in &pedigree.individuals {
        let mut marker_map = std::collections::HashMap::new();
        for marker in &marker_set.markers {
            let alleles = individual.get_alleles_by_marker_name(&marker.name);
            if !alleles.is_empty() {
                let allele_str = alleles
                    .iter()
                    .map(|a| a.display())
                    .collect::<Vec<_>>()
                    .join(";");
                marker_map.insert(marker.name.clone(), allele_str);
            }
        }
        if !marker_map.is_empty() {
            table.insert(individual.name.clone(), marker_map);
        }
    }

    let trace_alleles = trace.map(|h| {
        let mut marker_map = std::collections::HashMap::new();
        for marker in &marker_set.markers {
            let alleles = h.get_alleles_by_marker_name(&marker.name);
            if !alleles.is_empty() {
                let allele_str = alleles
                    .iter()
                    .map(|a| a.display())
                    .collect::<Vec<_>>()
                    .join(";");
                marker_map.insert(marker.name.clone(), allele_str);
            }
        }
        marker_map
    });

    Ok(HaplotypesParseResult {
        haplotype_table: table,
        trace_haplotype: trace_alleles,
        marker_names: marker_set.markers.iter().map(|m| m.name.clone()).collect(),
    })
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct HaplotypesParseResult {
    /// individual_name → marker_name → allele_string
    pub haplotype_table:
        std::collections::HashMap<String, std::collections::HashMap<String, String>>,
    /// marker_name → allele_string (for TRACE profile)
    pub trace_haplotype: Option<std::collections::HashMap<String, String>>,
    /// Ordered list of marker names
    pub marker_names: Vec<String>,
}

/// Serialize the haplotype table back to a JSON string.
#[tauri::command]
pub fn export_haplotypes_json(
    haplotype_table: std::collections::HashMap<String, std::collections::HashMap<String, String>>,
    trace_haplotype: Option<std::collections::HashMap<String, String>>,
) -> Result<String, String> {
    let mut out: serde_json::Map<String, Value> = serde_json::Map::new();
    if let Some(trace) = trace_haplotype {
        out.insert(
            "TRACE".into(),
            serde_json::to_value(trace).map_err(|e| e.to_string())?,
        );
    }
    for (name, markers) in haplotype_table {
        out.insert(
            name,
            serde_json::to_value(markers).map_err(|e| e.to_string())?,
        );
    }
    serde_json::to_string_pretty(&Value::Object(out)).map_err(|e| e.to_string())
}
