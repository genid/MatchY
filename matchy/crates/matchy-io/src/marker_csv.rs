/// Marker set CSV reader and writer.
///
/// Format:
///   marker,mutation_rate
///   DYS576,0.0021
///   DYS390,0.00001
use crate::{IoError, Result};
use matchy_core::{Marker, MarkerSet};
use std::io::{Read, Write};

/// Parse a CSV marker set.
pub fn read_marker_csv<R: Read>(reader: R) -> Result<MarkerSet> {
    let mut rdr = csv::Reader::from_reader(reader);
    let mut marker_set = MarkerSet::new();

    for result in rdr.records() {
        let record = result?;
        if record.len() < 2 {
            return Err(IoError::InvalidFormat(format!(
                "Expected 2 columns in marker CSV, got {}",
                record.len()
            )));
        }
        let name = record[0].trim();
        let rate_str = record[1].trim();

        if name.is_empty() {
            tracing::warn!("Empty marker name in CSV — skipping");
            continue;
        }

        let mutation_rate = rate_str.parse::<f64>().map_err(|_| {
            IoError::InvalidFormat(format!("Invalid mutation rate '{}' for marker '{}'", rate_str, name))
        })?;

        if mutation_rate < 0.0 || mutation_rate > 1.0 {
            return Err(IoError::InvalidFormat(format!(
                "Mutation rate {} out of [0, 1] for marker '{}'",
                mutation_rate, name
            )));
        }

        marker_set.add_marker(Marker::new(name, mutation_rate));
    }

    Ok(marker_set)
}

/// Write a MarkerSet to CSV format.
pub fn write_marker_csv<W: Write>(marker_set: &MarkerSet, writer: W) -> Result<()> {
    let mut wtr = csv::Writer::from_writer(writer);
    wtr.write_record(["marker", "mutation_rate"])?;
    for marker in &marker_set.markers {
        wtr.write_record([&marker.name, &marker.mutation_rate.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convenience: read from a string.
pub fn read_marker_csv_str(s: &str) -> Result<MarkerSet> {
    read_marker_csv(s.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_csv() {
        let input = "marker,mutation_rate\nDYS576,0.0021\nDYS390,0.00001\n";
        let ms = read_marker_csv_str(input).unwrap();
        assert_eq!(ms.markers.len(), 2);
        assert_eq!(ms.markers[0].name, "DYS576");
        assert!((ms.markers[0].mutation_rate - 0.0021).abs() < 1e-10);
    }
}
