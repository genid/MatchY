/// Haplotype JSON reader and writer.
///
/// Format:
/// {
///   "TRACE": { "DYS576": "15;16", ... },   // optional trace profile
///   "IndividualName": { "DYS576": "15", ... }
/// }
///
/// Allele strings: "15", "15.2" (with intermediate), "15;16" (multi-copy)
use crate::{IoError, Result};
use matchy_core::{Haplotype, HaplotypeClass, MarkerSet, Pedigree};
use serde_json::Value;
use std::collections::HashMap;
use std::io::Read;

/// Parse allele string into (value, intermediate_value) pairs.
/// "15"   → [(15, None)]
/// "15.2" → [(15, Some(2))]
/// "15;16.3" → [(15, None), (16, Some(3))]
fn parse_allele_string(s: &str) -> Result<Vec<(i32, Option<i32>)>> {
    s.split(';')
        .map(|part| {
            let part = part.trim();
            if let Some(dot_pos) = part.find('.') {
                let value_str = &part[..dot_pos];
                let iv_str = &part[dot_pos + 1..];
                let value = value_str.parse::<i32>().map_err(|_| {
                    IoError::InvalidFormat(format!("Invalid allele value: '{}'", value_str))
                })?;
                let iv = iv_str.parse::<i32>().map_err(|_| {
                    IoError::InvalidFormat(format!("Invalid intermediate value: '{}'", iv_str))
                })?;
                Ok((value, Some(iv)))
            } else {
                let value = part.parse::<i32>().map_err(|_| {
                    IoError::InvalidFormat(format!("Invalid allele value: '{}'", part))
                })?;
                Ok((value, None))
            }
        })
        .collect()
}

/// Serialize allele (value, intermediate_value) pairs back to string.
fn format_allele_string(alleles: &[(i32, Option<i32>)]) -> String {
    alleles
        .iter()
        .map(|(v, iv)| match iv {
            None => v.to_string(),
            Some(i) => format!("{}.{}", v, i),
        })
        .collect::<Vec<_>>()
        .join(";")
}

/// Read haplotypes from a JSON reader and assign them to individuals in the pedigree.
/// Returns the trace haplotype if a "TRACE" key was present, otherwise None.
pub fn read_haplotypes<R: Read>(
    reader: R,
    pedigree: &mut Pedigree,
    marker_set: &mut MarkerSet,
) -> Result<Option<Haplotype>> {
    let json: serde_json::Value = serde_json::from_reader(reader)?;
    let obj = json
        .as_object()
        .ok_or_else(|| IoError::InvalidFormat("Expected JSON object at top level".into()))?;

    let mut trace_haplotype: Option<Haplotype> = None;

    for (individual_name, marker_alleles_val) in obj {
        let marker_alleles = marker_alleles_val.as_object().ok_or_else(|| {
            IoError::InvalidFormat(format!(
                "Expected object for individual '{}'",
                individual_name
            ))
        })?;
        let mut haplotype = Haplotype::new();

        for (marker_name, allele_val) in marker_alleles {
            let allele_string = allele_val.as_str().ok_or_else(|| {
                IoError::InvalidFormat(format!(
                    "Expected string allele for marker '{}' in individual '{}'",
                    marker_name, individual_name
                ))
            })?;
            let allele_pairs = parse_allele_string(allele_string)?;
            let copy_count = allele_pairs.len() as u32;

            // Ensure marker exists in set; update copy count if needed
            if let Some(marker) = marker_set.markers.iter_mut().find(|m| m.name == *marker_name) {
                if marker.number_of_copies.is_none() {
                    marker.number_of_copies = Some(copy_count);
                } else if marker.number_of_copies != Some(copy_count) {
                    tracing::error!(
                        "Copy count mismatch for marker {} in individual {}",
                        marker_name,
                        individual_name
                    );
                    continue;
                }
                let marker_clone = marker.clone();
                for (value, iv) in &allele_pairs {
                    haplotype.add_allele(marker_clone.clone(), *value, *iv);
                }
            } else {
                tracing::warn!(
                    "Marker '{}' not found in marker set — skipping",
                    marker_name
                );
                continue;
            }
        }

        if individual_name == "TRACE" {
            // Check no pedigree individual is named TRACE
            if pedigree.get_individual_by_name("TRACE").is_some() {
                tracing::error!(
                    "Ambiguity: JSON contains 'TRACE' key but pedigree has individual named 'TRACE'"
                );
                return Err(IoError::InvalidFormat(
                    "Individual named 'TRACE' conflicts with reserved TRACE key".into(),
                ));
            }
            trace_haplotype = Some(haplotype.clone());
            pedigree.trace_haplotype = Some(haplotype);
        } else if let Some(individual) = pedigree.get_individual_by_name_mut(individual_name) {
            individual.haplotype = haplotype;
            individual.haplotype_class = HaplotypeClass::Known;
        } else {
            tracing::warn!(
                "Individual '{}' not found in pedigree — skipping haplotype",
                individual_name
            );
        }
    }

    Ok(trace_haplotype)
}

/// Serialize pedigree haplotypes to a JSON value.
pub fn write_haplotypes(
    pedigree: &Pedigree,
    trace_haplotype: Option<&Haplotype>,
    marker_set: &MarkerSet,
) -> Value {
    let mut out: HashMap<String, HashMap<String, String>> = HashMap::new();

    if let Some(trace) = trace_haplotype {
        let mut marker_map = HashMap::new();
        for marker in &marker_set.markers {
            let alleles: Vec<(i32, Option<i32>)> = trace
                .get_alleles_by_marker_name(&marker.name)
                .iter()
                .map(|a| (a.value, a.intermediate_value))
                .collect();
            if !alleles.is_empty() {
                marker_map.insert(marker.name.clone(), format_allele_string(&alleles));
            }
        }
        out.insert("TRACE".into(), marker_map);
    }

    for individual in &pedigree.individuals {
        if individual.haplotype.alleles.is_empty() {
            continue;
        }
        let mut marker_map = HashMap::new();
        for marker in &marker_set.markers {
            let alleles: Vec<(i32, Option<i32>)> = individual
                .get_alleles_by_marker_name(&marker.name)
                .iter()
                .map(|a| (a.value, a.intermediate_value))
                .collect();
            if !alleles.is_empty() {
                marker_map.insert(marker.name.clone(), format_allele_string(&alleles));
            }
        }
        if !marker_map.is_empty() {
            out.insert(individual.name.clone(), marker_map);
        }
    }

    serde_json::to_value(out).unwrap_or(Value::Null)
}

#[cfg(test)]
mod tests {
    use super::*;
    use matchy_core::{Marker, MarkerSet, Pedigree};

    fn make_marker_set() -> MarkerSet {
        let mut ms = MarkerSet::new();
        ms.add_marker(Marker::new("DYS576", 0.0021));
        ms.add_marker(Marker::new("DYS390", 0.00001));
        ms
    }

    #[test]
    fn parse_simple_allele() {
        let pairs = parse_allele_string("15").unwrap();
        assert_eq!(pairs, vec![(15, None)]);
    }

    #[test]
    fn parse_intermediate_allele() {
        let pairs = parse_allele_string("15.2").unwrap();
        assert_eq!(pairs, vec![(15, Some(2))]);
    }

    #[test]
    fn parse_multicopy_allele() {
        let pairs = parse_allele_string("15;16").unwrap();
        assert_eq!(pairs, vec![(15, None), (16, None)]);
    }

    #[test]
    fn read_basic_haplotypes() {
        let json_str = r#"{"suspect": {"DYS576": "18", "DYS390": "24"}}"#;
        let mut pedigree = Pedigree::new();
        pedigree.add_individual("1", "suspect");
        let mut ms = make_marker_set();
        let trace = read_haplotypes(json_str.as_bytes(), &mut pedigree, &mut ms).unwrap();
        assert!(trace.is_none());
        let ind = pedigree.get_individual_by_name("suspect").unwrap();
        assert_eq!(ind.haplotype_class, matchy_core::HaplotypeClass::Known);
    }
}
