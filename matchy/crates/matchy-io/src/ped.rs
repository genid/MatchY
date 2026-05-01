/// PLINK/LINKAGE .ped format reader.
///
/// Columns: family_id individual_id paternal_id maternal_id sex phenotype
/// Only male individuals (sex == 1) are included.
/// Females (sex == 2) are skipped.
use crate::{IoError, Result};
use matchy_core::Pedigree;
use std::io::BufRead;

/// Parse a .ped file into a Pedigree (males only).
pub fn read_ped<R: BufRead>(reader: R) -> Result<Pedigree> {
    let mut pedigree = Pedigree::new();
    let mut paternal_links: Vec<(String, String)> = Vec::new(); // (paternal_id, child_id)

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 6 {
            return Err(IoError::InvalidFormat(format!(
                "Expected 6 columns in .ped line, got {}: '{}'",
                parts.len(),
                line
            )));
        }
        let (_family_id, individual_id, paternal_id, _maternal_id, sex, _phenotype) =
            (parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]);

        if sex == "2" {
            continue; // Skip females
        }
        if sex != "1" {
            tracing::warn!("Unexpected sex code '{}' in .ped file — skipping", sex);
            continue;
        }

        pedigree.add_individual(individual_id, individual_id);

        if paternal_id != "0" {
            paternal_links.push((paternal_id.to_string(), individual_id.to_string()));
        }
    }

    for (parent_id, child_id) in paternal_links {
        // Only add relationship if both parent and child exist in the pedigree
        if pedigree.get_individual_by_id(&parent_id).is_some()
            && pedigree.get_individual_by_id(&child_id).is_some()
        {
            pedigree.add_relationship(parent_id, child_id);
        }
    }

    Ok(pedigree)
}

/// Convenience: read from a string.
pub fn read_ped_str(s: &str) -> Result<Pedigree> {
    read_ped(std::io::BufReader::new(s.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_ped() {
        let input = "FAM1 1 0 0 1 0\nFAM1 2 1 0 1 0\nFAM1 3 0 0 2 0\n";
        let pedigree = read_ped_str(input).unwrap();
        // Individual 3 is female — excluded
        assert_eq!(pedigree.individuals.len(), 2);
        assert_eq!(pedigree.relationships.len(), 1);
        assert_eq!(pedigree.relationships[0].parent_id, "1");
        assert_eq!(pedigree.relationships[0].child_id, "2");
    }
}
