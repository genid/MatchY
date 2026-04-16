/// Trivial Graph Format (.tgf) reader and writer.
///
/// Format:
///   <id> <name>        (one line per individual, name may equal id)
///   #                  (separator)
///   <parent_id> <child_id>  (one line per relationship)
use crate::{IoError, Result};
use matchy_core::Pedigree;
use std::io::{BufRead, Write};

/// Parse a .tgf file into a Pedigree.
pub fn read_tgf<R: BufRead>(reader: R) -> Result<Pedigree> {
    let mut pedigree = Pedigree::new();
    let mut in_edges = false;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line == "#" {
            in_edges = true;
            continue;
        }
        if in_edges {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 2 {
                return Err(IoError::InvalidFormat(format!(
                    "Invalid edge line in TGF: '{}'",
                    line
                )));
            }
            pedigree.add_relationship(parts[0], parts[1]);
        } else {
            let parts: Vec<&str> = line.splitn(2, ' ').collect();
            match parts.len() {
                2 => pedigree.add_individual(parts[0], parts[1]),
                1 => pedigree.add_individual(parts[0], parts[0]),
                _ => {
                    return Err(IoError::InvalidFormat(format!(
                        "Invalid node line in TGF: '{}'",
                        line
                    )))
                }
            }
        }
    }
    Ok(pedigree)
}

/// Serialize a Pedigree to TGF format.
pub fn write_tgf<W: Write>(pedigree: &Pedigree, writer: &mut W) -> Result<()> {
    for individual in &pedigree.individuals {
        writeln!(writer, "{} {}", individual.id, individual.name)?;
    }
    writeln!(writer, "#")?;
    for rel in &pedigree.relationships {
        writeln!(writer, "{} {}", rel.parent_id, rel.child_id)?;
    }
    Ok(())
}

/// Convenience: read TGF from a string.
pub fn read_tgf_str(s: &str) -> Result<Pedigree> {
    read_tgf(std::io::BufReader::new(s.as_bytes()))
}

/// Convenience: write TGF to a String.
pub fn write_tgf_str(pedigree: &Pedigree) -> Result<String> {
    let mut buf = Vec::new();
    write_tgf(pedigree, &mut buf)?;
    Ok(String::from_utf8(buf).unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_tgf() {
        let input = "1 grandfather\n2 father\n3 suspect\n#\n1 2\n2 3\n";
        let pedigree = read_tgf_str(input).unwrap();
        assert_eq!(pedigree.individuals.len(), 3);
        assert_eq!(pedigree.relationships.len(), 2);

        let output = write_tgf_str(&pedigree).unwrap();
        let pedigree2 = read_tgf_str(&output).unwrap();
        assert_eq!(pedigree2.individuals.len(), 3);
        assert_eq!(pedigree2.relationships.len(), 2);
    }

    #[test]
    fn tgf_without_names() {
        let input = "1\n2\n#\n1 2\n";
        let pedigree = read_tgf_str(input).unwrap();
        assert_eq!(pedigree.individuals[0].name, "1");
    }
}
