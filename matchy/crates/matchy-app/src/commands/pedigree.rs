use matchy_core::Pedigree;
use serde::{Deserialize, Serialize};

/// Serializable pedigree representation for the frontend.
#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PedigreeData {
    pub individuals: Vec<IndividualData>,
    pub relationships: Vec<RelationshipData>,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IndividualData {
    pub id: String,
    pub name: String,
    pub haplotype_class: String,
    pub exclude: bool,
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RelationshipData {
    pub parent_id: String,
    pub child_id: String,
}

fn pedigree_to_data(pedigree: &Pedigree) -> PedigreeData {
    PedigreeData {
        individuals: pedigree
            .individuals
            .iter()
            .map(|i| IndividualData {
                id: i.id.clone(),
                name: i.name.clone(),
                haplotype_class: i.haplotype_class.as_str().to_string(),
                exclude: i.exclude,
            })
            .collect(),
        relationships: pedigree
            .relationships
            .iter()
            .map(|r| RelationshipData {
                parent_id: r.parent_id.clone(),
                child_id: r.child_id.clone(),
            })
            .collect(),
    }
}

/// Parse a TGF string into a pedigree structure.
#[tauri::command]
pub fn parse_tgf(tgf_content: String) -> Result<PedigreeData, String> {
    let pedigree = matchy_io::tgf::read_tgf_str(&tgf_content).map_err(|e| e.to_string())?;
    Ok(pedigree_to_data(&pedigree))
}

/// Parse a PED string into a pedigree structure.
#[tauri::command]
pub fn parse_ped(ped_content: String) -> Result<PedigreeData, String> {
    let pedigree = matchy_io::ped::read_ped_str(&ped_content).map_err(|e| e.to_string())?;
    Ok(pedigree_to_data(&pedigree))
}

/// Serialize a pedigree structure to TGF string.
#[tauri::command]
pub fn export_tgf(data: PedigreeData) -> Result<String, String> {
    let mut pedigree = Pedigree::new();
    for ind in &data.individuals {
        pedigree.add_individual(ind.id.clone(), ind.name.clone());
    }
    for rel in &data.relationships {
        pedigree.add_relationship(rel.parent_id.clone(), rel.child_id.clone());
    }
    matchy_io::tgf::write_tgf_str(&pedigree).map_err(|e| e.to_string())
}

/// Validate a pedigree (DAG check, connectivity).
#[tauri::command]
pub fn validate_pedigree(data: PedigreeData) -> Result<ValidationResult, String> {
    let mut pedigree = Pedigree::new();
    for ind in &data.individuals {
        pedigree.add_individual(ind.id.clone(), ind.name.clone());
    }
    for rel in &data.relationships {
        pedigree.add_relationship(rel.parent_id.clone(), rel.child_id.clone());
    }

    let dag_ok = matchy_core::graph::validate_dag(&pedigree).is_ok();
    let connected = matchy_core::graph::is_connected(&pedigree);

    Ok(ValidationResult {
        is_valid: dag_ok && connected,
        is_dag: dag_ok,
        is_connected: connected,
        error: if !dag_ok {
            Some("Pedigree contains cycles".into())
        } else if !connected {
            Some("Pedigree is not connected".into())
        } else {
            None
        },
    })
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ValidationResult {
    pub is_valid: bool,
    pub is_dag: bool,
    pub is_connected: bool,
    pub error: Option<String>,
}
