/// Embedded mutation rates database.
///
/// Loaded from the bundled mutation_rates.csv (embedded at compile time).
/// Used to look up rates for markers not in a custom file.
use crate::Result;
use std::collections::HashMap;

const MUTATION_RATES_CSV: &str = include_str!("../assets/mutation_rates.csv");

/// Load all mutation rates from the embedded CSV.
/// Returns a map: marker_name → mutation_rate.
pub fn load_mutation_rates() -> Result<HashMap<String, f64>> {
    let mut rdr = csv::Reader::from_reader(MUTATION_RATES_CSV.as_bytes());
    let mut rates = HashMap::new();
    for record in rdr.records() {
        let record = record?;
        if record.len() < 2 {
            continue;
        }
        let name = record[0].trim().to_string();
        let rate_str = record[1].trim();
        if let Ok(rate) = rate_str.parse::<f64>() {
            rates.insert(name, rate);
        }
    }
    Ok(rates)
}

/// Look up mutation rate for a single marker.
pub fn get_mutation_rate(marker_name: &str) -> Result<Option<f64>> {
    let rates = load_mutation_rates()?;
    Ok(rates.get(marker_name).copied())
}
