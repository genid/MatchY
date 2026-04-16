use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use rust_decimal::Decimal;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Marker
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Marker {
    pub name: String,
    pub mutation_rate: f64,
    /// Number of copies (multi-copy markers like DYF387S1 have 2).
    /// None means single-copy (treated as 1).
    pub number_of_copies: Option<u32>,
}

impl Marker {
    pub fn new(name: impl Into<String>, mutation_rate: f64) -> Self {
        Self {
            name: name.into(),
            mutation_rate,
            number_of_copies: None,
        }
    }

    pub fn with_copies(mut self, n: u32) -> Self {
        self.number_of_copies = Some(n);
        self
    }

    pub fn copies(&self) -> u32 {
        self.number_of_copies.unwrap_or(1)
    }

    /// Single-copy mutation rate derived from the all-copy rate:
    ///   mu_single = 1 - (1 - mu_all)^(1/n)
    pub fn single_copy_mutation_rate(&self) -> f64 {
        let n = self.copies() as f64;
        if n <= 1.0 {
            self.mutation_rate
        } else {
            1.0 - (1.0 - self.mutation_rate).powf(1.0 / n)
        }
    }
}

impl PartialEq for Marker {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
            && self.mutation_rate == other.mutation_rate
            && self.number_of_copies == other.number_of_copies
    }
}

impl Eq for Marker {}

impl Hash for Marker {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name.hash(state);
        // Convert f64 to bits for hashing (consistent with Python's float hash for finite values)
        self.mutation_rate.to_bits().hash(state);
        self.number_of_copies.hash(state);
    }
}

// ---------------------------------------------------------------------------
// Allele
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Allele {
    pub marker: Marker,
    pub value: i32,
    /// Fractional sub-unit (e.g. the "2" in "15.2"). None = no intermediate.
    pub intermediate_value: Option<i32>,
    /// Number of mutation steps from parent allele (set during simulation).
    pub mutation_value: Option<i32>,
    /// Probability of mutation (set during simulation).
    pub mutation_probability: Option<f64>,
}

impl Allele {
    pub fn new(marker: Marker, value: i32) -> Self {
        Self {
            marker,
            value,
            intermediate_value: None,
            mutation_value: None,
            mutation_probability: None,
        }
    }

    pub fn with_intermediate(mut self, iv: i32) -> Self {
        self.intermediate_value = Some(iv);
        self
    }

    /// Normalised intermediate value: treats None as 0, matching Python's hash.
    pub fn intermediate_norm(&self) -> i32 {
        self.intermediate_value.unwrap_or(0)
    }

    /// Display string: "15" or "15.2"
    pub fn display(&self) -> String {
        match self.intermediate_value {
            None => self.value.to_string(),
            Some(iv) => format!("{}.{}", self.value, iv),
        }
    }
}

impl PartialEq for Allele {
    fn eq(&self, other: &Self) -> bool {
        self.marker == other.marker
            && self.value == other.value
            && self.intermediate_value == other.intermediate_value
    }
}

impl Eq for Allele {}

impl Hash for Allele {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.marker.hash(state);
        self.value.hash(state);
        self.intermediate_value.hash(state);
    }
}

impl PartialOrd for Allele {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Allele {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (self.value, self.intermediate_norm()).cmp(&(other.value, other.intermediate_norm()))
    }
}

// ---------------------------------------------------------------------------
// Haplotype
// ---------------------------------------------------------------------------

/// A collection of alleles grouped by marker name.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Haplotype {
    /// marker_name → sorted list of alleles (multi-copy markers have >1 allele)
    pub alleles: HashMap<String, Vec<Allele>>,
}

impl Haplotype {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_allele(&mut self, marker: Marker, value: i32, intermediate_value: Option<i32>) {
        let allele = if let Some(iv) = intermediate_value {
            Allele::new(marker.clone(), value).with_intermediate(iv)
        } else {
            Allele::new(marker.clone(), value)
        };
        self.alleles
            .entry(marker.name.clone())
            .or_default()
            .push(allele);
    }

    pub fn get_alleles_by_marker_name(&self, marker_name: &str) -> Vec<&Allele> {
        let mut v: Vec<&Allele> = self.alleles.get(marker_name).map(|a| a.iter().collect()).unwrap_or_default();
        v.sort();
        v
    }

    /// Allelic distance to another haplotype (sum of abs value differences, paired by sorted order).
    pub fn allelic_difference(&self, other: &Haplotype) -> Option<i32> {
        if self.alleles.keys().collect::<std::collections::HashSet<_>>()
            != other.alleles.keys().collect::<std::collections::HashSet<_>>()
        {
            return None;
        }
        let mut diff = 0i32;
        for (marker_name, alleles) in &self.alleles {
            let other_alleles = other.alleles.get(marker_name)?;
            let mut a: Vec<&Allele> = alleles.iter().collect();
            let mut b: Vec<&Allele> = other_alleles.iter().collect();
            a.sort();
            b.sort();
            for (x, y) in a.iter().zip(b.iter()) {
                diff += (x.value - y.value).abs();
            }
        }
        Some(diff)
    }
}

impl PartialEq for Haplotype {
    fn eq(&self, other: &Self) -> bool {
        if self.alleles.len() != other.alleles.len() {
            return false;
        }
        for (marker, alleles) in &self.alleles {
            let Some(other_alleles) = other.alleles.get(marker) else {
                return false;
            };
            if alleles.len() != other_alleles.len() {
                return false;
            }
            let mut a: Vec<&Allele> = alleles.iter().collect();
            let mut b: Vec<&Allele> = other_alleles.iter().collect();
            a.sort();
            b.sort();
            for (x, y) in a.iter().zip(b.iter()) {
                if x.value != y.value || x.intermediate_value != y.intermediate_value {
                    return false;
                }
            }
        }
        true
    }
}

impl Eq for Haplotype {}

impl Hash for Haplotype {
    /// Deterministic hash matching Python's Haplotype.__hash__:
    /// Sort by marker name, sort alleles within each marker,
    /// normalise intermediate_value=None → 0.
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut sorted_markers: Vec<(&String, &Vec<Allele>)> = self.alleles.iter().collect();
        sorted_markers.sort_by_key(|(k, _)| k.as_str());
        for (marker_name, alleles) in sorted_markers {
            marker_name.hash(state);
            let mut sorted_alleles: Vec<&Allele> = alleles.iter().collect();
            sorted_alleles.sort();
            for allele in sorted_alleles {
                allele.value.hash(state);
                allele.intermediate_norm().hash(state);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// HaplotypeClass
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
#[serde(rename_all = "lowercase")]
pub enum HaplotypeClass {
    #[default]
    Unknown,
    Known,
    Suspect,
    Estimated,
    Fixed,
    Excluded,
}

impl HaplotypeClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            HaplotypeClass::Unknown => "unknown",
            HaplotypeClass::Known => "known",
            HaplotypeClass::Suspect => "suspect",
            HaplotypeClass::Estimated => "estimated",
            HaplotypeClass::Fixed => "fixed",
            HaplotypeClass::Excluded => "excluded",
        }
    }
}

impl std::fmt::Display for HaplotypeClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

// ---------------------------------------------------------------------------
// Individual
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    pub id: String,
    pub name: String,
    pub haplotype: Haplotype,
    pub haplotype_class: HaplotypeClass,
    pub exclude: bool,
    pub picking_probability: Option<Decimal>,
    pub closest_known_individuals: Vec<String>, // store IDs
    pub closest_known_distance: Option<u32>,
}

impl Individual {
    pub fn new(id: impl Into<String>, name: impl Into<String>) -> Self {
        Self {
            id: id.into(),
            name: name.into(),
            haplotype: Haplotype::new(),
            haplotype_class: HaplotypeClass::Unknown,
            exclude: false,
            picking_probability: None,
            closest_known_individuals: Vec::new(),
            closest_known_distance: None,
        }
    }

    pub fn add_allele(&mut self, marker: Marker, value: i32, intermediate_value: Option<i32>) {
        self.haplotype.add_allele(marker, value, intermediate_value);
    }

    pub fn get_alleles_by_marker_name(&self, marker_name: &str) -> Vec<&Allele> {
        self.haplotype.get_alleles_by_marker_name(marker_name)
    }

    pub fn is_known(&self) -> bool {
        matches!(
            self.haplotype_class,
            HaplotypeClass::Known | HaplotypeClass::Suspect | HaplotypeClass::Fixed
        )
    }
}

impl PartialEq for Individual {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Individual {}

impl Hash for Individual {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

// ---------------------------------------------------------------------------
// Relationship
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Relationship {
    pub parent_id: String,
    pub child_id: String,
    pub edge_class: String,
}

impl Relationship {
    pub fn new(parent_id: impl Into<String>, child_id: impl Into<String>) -> Self {
        Self {
            parent_id: parent_id.into(),
            child_id: child_id.into(),
            edge_class: "unknown".into(),
        }
    }
}

// ---------------------------------------------------------------------------
// MarkerSet
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct MarkerSet {
    pub markers: Vec<Marker>,
}

impl MarkerSet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_marker(&mut self, marker: Marker) {
        self.markers.push(marker);
    }

    pub fn get_marker_by_name(&self, name: &str) -> Option<&Marker> {
        self.markers.iter().find(|m| m.name == name)
    }
}

impl Hash for MarkerSet {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut sorted: Vec<_> = self.markers.iter().map(|m| (&m.name, &m.mutation_rate, &m.number_of_copies)).collect();
        sorted.sort_by_key(|(name, _, _)| name.as_str());
        for (name, rate, copies) in sorted {
            name.hash(state);
            rate.to_bits().hash(state);
            copies.hash(state);
        }
    }
}

// ---------------------------------------------------------------------------
// Pedigree
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Pedigree {
    pub individuals: Vec<Individual>,
    pub relationships: Vec<Relationship>,
    pub picking_probabilities: HashMap<String, Decimal>,
}

impl Pedigree {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_individual(&mut self, id: impl Into<String>, name: impl Into<String>) {
        let id = id.into();
        let name = name.into();
        if self.individuals.iter().any(|i| i.name == name) {
            tracing::warn!("Individual with name '{}' already exists", name);
        }
        self.individuals.push(Individual::new(id, name));
    }

    pub fn add_relationship(&mut self, parent_id: impl Into<String>, child_id: impl Into<String>) {
        self.relationships.push(Relationship::new(parent_id, child_id));
    }

    pub fn get_individual_by_id(&self, id: &str) -> Option<&Individual> {
        self.individuals.iter().find(|i| i.id == id)
    }

    pub fn get_individual_by_id_mut(&mut self, id: &str) -> Option<&mut Individual> {
        self.individuals.iter_mut().find(|i| i.id == id)
    }

    pub fn get_individual_by_name(&self, name: &str) -> Option<&Individual> {
        self.individuals.iter().find(|i| i.name == name)
    }

    pub fn get_individual_by_name_mut(&mut self, name: &str) -> Option<&mut Individual> {
        self.individuals.iter_mut().find(|i| i.name == name)
    }

    /// Returns parent ID strings for a given individual ID.
    pub fn parents_of(&self, child_id: &str) -> Vec<&str> {
        self.relationships
            .iter()
            .filter(|r| r.child_id == child_id)
            .map(|r| r.parent_id.as_str())
            .collect()
    }

    /// Returns child ID strings for a given individual ID.
    pub fn children_of(&self, parent_id: &str) -> Vec<&str> {
        self.relationships
            .iter()
            .filter(|r| r.parent_id == parent_id)
            .map(|r| r.child_id.as_str())
            .collect()
    }

    /// Returns IDs of all individuals with no parents (roots of the DAG).
    pub fn roots(&self) -> Vec<&str> {
        self.individuals
            .iter()
            .filter(|i| self.parents_of(&i.id).is_empty())
            .map(|i| i.id.as_str())
            .collect()
    }
}

impl Hash for Pedigree {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut ind_tuples: Vec<_> = self
            .individuals
            .iter()
            .map(|i| (i.id.as_str(), i.name.as_str(), i.haplotype_class.as_str()))
            .collect();
        ind_tuples.sort();
        ind_tuples.hash(state);

        let mut rel_tuples: Vec<_> = self
            .relationships
            .iter()
            .map(|r| (r.parent_id.as_str(), r.child_id.as_str()))
            .collect();
        rel_tuples.sort();
        rel_tuples.hash(state);
    }
}

// ---------------------------------------------------------------------------
// Bias
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bias {
    pub marker: Marker,
    pub copy_nr: u32,
    pub direction: BiasDirection,
    pub target_mass: f64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum BiasDirection {
    Up,
    Down,
}

// ---------------------------------------------------------------------------
// SimulationParameters
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationParameters {
    pub two_step_mutation_fraction: f64,
    pub batch_length: u64,
    pub convergence_criterion: f64,
    /// None = auto-calculate per-individual
    pub bias: Option<f64>,
    pub number_of_threads: usize,
    /// Suspect individual name (None = trace mode)
    pub suspect: Option<String>,
    pub exclude: Vec<String>,
    pub skip_inside: bool,
    pub skip_outside: bool,
    pub trace_mode: bool,
    pub adaptive_bias: bool,
    pub simulation_name: String,
    pub user_name: String,
    pub results_path: std::path::PathBuf,
}

impl Default for SimulationParameters {
    fn default() -> Self {
        Self {
            two_step_mutation_fraction: 0.03,
            batch_length: 10_000,
            convergence_criterion: 0.02,
            bias: None,
            number_of_threads: 4,
            suspect: None,
            exclude: Vec::new(),
            skip_inside: false,
            skip_outside: false,
            trace_mode: false,
            adaptive_bias: false,
            simulation_name: "simulation".into(),
            user_name: String::new(),
            results_path: std::path::PathBuf::from("./results"),
        }
    }
}

// ---------------------------------------------------------------------------
// SimulationResult
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    pub parameters: SimulationParameters,
    /// P(observed | inside match probability at x matches)
    pub inside_match_probabilities: Option<MatchProbabilities>,
    /// P(random outside male matches suspect)
    pub outside_match_probability: Option<Decimal>,
    /// Per-individual marginal probabilities (trace mode)
    pub per_individual_probabilities: Option<HashMap<String, Decimal>>,
    pub trials: u32,
    pub converged: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchProbabilities {
    /// match_count → probability
    pub probabilities: HashMap<u32, Decimal>,
    pub average_pedigree_probability: Decimal,
}
