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
    /// Crime scene trace haplotype (loaded from the "TRACE" key in the JSON).
    /// Used as the comparison profile in trace mode.
    #[serde(skip)]
    pub trace_haplotype: Option<Haplotype>,
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

    /// Extend the pedigree for the outside-match calculation.
    ///
    /// Mirrors Python's `extend_pedigree()` in models.py:700–730.
    ///
    /// Adds one unknown parent above the current root, then a chain of
    /// `highest_level_with_known` unknown children below it.
    /// `highest_level_with_known` is the 0-indexed BFS depth of the deepest
    /// known individual in the current pedigree.
    ///
    /// Returns the ID of the last added child — the "outside" individual.
    pub fn extend_pedigree(&mut self) -> String {
        use crate::graph::bfs_layers;

        // Find the current root.
        let root_ids = self.roots();
        let root_id = root_ids
            .first()
            .expect("extend_pedigree: pedigree has no root")
            .to_string();

        // BFS layers from the root.
        // Mirror Python models.py:700-730: find the FIRST (shallowest) layer that
        // contains a "known" individual (not "suspect" or other classes), then set
        // highest = that_depth + 1.  If no known individual is found, use the total
        // number of layers.
        let layers = bfs_layers(self, &root_id).unwrap_or_default();
        let first_known_depth = layers
            .iter()
            .enumerate()
            .find_map(|(depth, layer)| {
                let has_known = layer.iter().any(|id| {
                    self.get_individual_by_id(id)
                        .map(|ind| ind.haplotype_class == HaplotypeClass::Known)
                        .unwrap_or(false)
                });
                if has_known { Some(depth) } else { None }
            });
        let highest_level_with_known = match first_known_depth {
            Some(d) => d + 1,
            None => layers.len(),
        };

        // Generate unique integer IDs not already in use.
        let existing_ids: std::collections::HashSet<u64> = self
            .individuals
            .iter()
            .filter_map(|i| i.id.parse::<u64>().ok())
            .collect();
        let mut next_id = (0u64..).find(|n| !existing_ids.contains(n)).unwrap_or(0);
        let mut alloc_id = || -> String {
            while existing_ids.contains(&next_id) {
                next_id += 1;
            }
            let id = next_id.to_string();
            next_id += 1;
            id
        };

        // Add new_root above old root.
        let new_root_id = alloc_id();
        self.add_individual(new_root_id.clone(), "new_root".to_string());
        self.add_relationship(new_root_id.clone(), root_id.clone());

        // Add a chain of unknowns below new_root.
        let mut prev_id = new_root_id.clone();
        let mut last_child_name = String::new();
        for i in 1..=highest_level_with_known {
            let child_id = alloc_id();
            let child_name = format!("new_child_{}", i);
            self.add_individual(child_id.clone(), child_name.clone());
            self.add_relationship(prev_id.clone(), child_id.clone());
            last_child_name = child_name;
            prev_id = child_id;
        }

        last_child_name
    }

    /// Reroot the pedigree DAG at `new_root_id`.
    ///
    /// Treats the pedigree as an undirected graph, performs a DFS from
    /// `new_root_id`, and replaces all relationships with the DFS spanning-tree
    /// edges directed away from the new root (parent → child).
    ///
    /// This mirrors Python's `reroot_pedigree`: after calling this, the named
    /// individual is the single root, and all edge directions follow the DFS tree.
    pub fn reroot(&mut self, new_root_id: &str) {
        use std::collections::{HashMap, HashSet};

        // Build undirected adjacency list
        let mut adj: HashMap<String, Vec<String>> = HashMap::new();
        for rel in &self.relationships {
            adj.entry(rel.parent_id.clone())
                .or_default()
                .push(rel.child_id.clone());
            adj.entry(rel.child_id.clone())
                .or_default()
                .push(rel.parent_id.clone());
        }

        // DFS from new_root_id — record parent_of[child] = parent
        let mut parent_of: HashMap<String, String> = HashMap::new();
        let mut visited: HashSet<String> = HashSet::new();
        let mut stack: Vec<String> = vec![new_root_id.to_string()];
        visited.insert(new_root_id.to_string());

        while let Some(node) = stack.pop() {
            if let Some(neighbours) = adj.get(&node) {
                for nbr in neighbours {
                    if !visited.contains(nbr.as_str()) {
                        visited.insert(nbr.clone());
                        parent_of.insert(nbr.clone(), node.clone());
                        stack.push(nbr.clone());
                    }
                }
            }
        }

        // Replace relationships with DFS tree edges (parent → child), sorted for determinism
        let mut rels: Vec<Relationship> = parent_of
            .into_iter()
            .map(|(child, parent)| Relationship::new(parent, child))
            .collect();
        rels.sort_by(|a, b| {
            a.parent_id
                .cmp(&b.parent_id)
                .then(a.child_id.cmp(&b.child_id))
        });
        self.relationships = rels;

        // Mirrors Python's reroot_pedigree: demote previous Suspect to Known,
        // then mark the new root as Suspect (used as anchor by calculate_picking_probabilities).
        for ind in &mut self.individuals {
            if ind.haplotype_class == HaplotypeClass::Suspect {
                ind.haplotype_class = HaplotypeClass::Known;
            }
        }
        if let Some(ind) = self.get_individual_by_id_mut(new_root_id) {
            if ind.haplotype_class == HaplotypeClass::Known {
                ind.haplotype_class = HaplotypeClass::Suspect;
            }
        }
    }

    /// Compute a priori picking probabilities for all unknown individuals.
    ///
    /// - Trace mode: uniform distribution (1/n for each unknown).
    /// - Suspect mode: mirrors Python models.py:653–697 exactly.
    ///   1. Compute shortest (undirected) paths from suspect to each unknown.
    ///   2. Estimate haplotypes for unknowns along those paths by copying from
    ///      their predecessor on the path.
    ///   3. For each unknown U: fix U to suspect haplotype, then sum allelic
    ///      differences of every unknown to its path-parent. (+1 to avoid zero.)
    ///   4. Normalise: (max - val + 1) / (max + 1).
    pub fn calculate_picking_probabilities(&mut self, trace_mode: bool) {
        let unknown_ids: Vec<String> = self
            .individuals
            .iter()
            .filter(|i| i.haplotype_class == HaplotypeClass::Unknown)
            .map(|i| i.id.clone())
            .collect();

        if unknown_ids.is_empty() {
            return;
        }

        if trace_mode {
            let prob = Decimal::ONE / Decimal::from(unknown_ids.len() as u32);
            self.picking_probabilities = unknown_ids
                .into_iter()
                .map(|id| (id, prob))
                .collect();
            return;
        }

        let suspect_id = match self.individuals.iter()
            .find(|i| i.haplotype_class == HaplotypeClass::Suspect)
        {
            Some(s) if !s.haplotype.alleles.is_empty() => s.id.clone(),
            _ => {
                tracing::error!("Suspect does not exist or has no haplotype.");
                return;
            }
        };
        let suspect_haplotype = self.get_individual_by_id(&suspect_id)
            .unwrap()
            .haplotype
            .clone();

        // Shortest paths from suspect to each unknown (undirected).
        let paths: HashMap<String, Vec<String>> = unknown_ids.iter()
            .filter_map(|uid| {
                let path = crate::graph::shortest_path_undirected(self, &suspect_id, uid)?;
                Some((uid.clone(), path))
            })
            .collect();

        // Clone pedigree and estimate haplotypes along paths (mirrors Python lines 668–675).
        let mut ped_est = self.clone();
        for uid in &unknown_ids {
            let path = match paths.get(uid) {
                Some(p) => p.clone(),
                None => continue,
            };
            for i in 0..path.len() {
                let node_id = &path[i];
                let class = ped_est.get_individual_by_id(node_id)
                    .map(|ind| ind.haplotype_class.clone());
                if class == Some(HaplotypeClass::Unknown) && i > 0 {
                    let prev_hap = ped_est.get_individual_by_id(&path[i - 1])
                        .map(|ind| ind.haplotype.clone());
                    if let Some(hap) = prev_hap {
                        if let Some(ind) = ped_est.get_individual_by_id_mut(node_id) {
                            ind.haplotype = hap;
                            ind.haplotype_class = HaplotypeClass::Estimated;
                        }
                    }
                }
            }
        }

        // For each unknown U, fix it to suspect and sum mutations across all unknowns.
        let mut scores: HashMap<String, i32> = HashMap::new();
        for uid in &unknown_ids {
            let mut ped_fixed = ped_est.clone();
            if let Some(ind) = ped_fixed.get_individual_by_id_mut(uid) {
                ind.haplotype = suspect_haplotype.clone();
                ind.haplotype_class = HaplotypeClass::Fixed;
            }
            let mut total: i32 = 0;
            for vid in &unknown_ids {
                let path = match paths.get(vid) {
                    Some(p) if p.len() >= 2 => p,
                    _ => continue,
                };
                let parent_id = &path[path.len() - 2];
                let v_hap = ped_fixed.get_individual_by_id(vid).map(|i| i.haplotype.clone());
                let p_hap = ped_fixed.get_individual_by_id(parent_id).map(|i| i.haplotype.clone());
                if let (Some(v), Some(p)) = (v_hap, p_hap) {
                    if let Some(diff) = v.allelic_difference(&p) {
                        total += diff;
                    }
                }
            }
            scores.insert(uid.clone(), total + 1);
        }

        // Normalise: (max - val + 1) / (max + 1).
        let max_val = scores.values().copied().max().unwrap_or(1);
        self.picking_probabilities = unknown_ids.iter()
            .filter_map(|uid| {
                let val = *scores.get(uid)?;
                let numerator = Decimal::from(max_val - val + 1);
                let denominator = Decimal::from(max_val + 1);
                Some((uid.clone(), numerator / denominator))
            })
            .collect();
    }

    /// Add a trace haplotype as a child of the first known individual.
    ///
    /// Mirrors Python's `add_trace()` in models.py:923-936.
    /// The new individual gets `haplotype_class = Suspect` so downstream logic
    /// (reroot, picking probs) treats it as the anchor.
    /// Returns the name of the newly added individual, or None if no known
    /// individual exists.
    pub fn add_trace(&mut self, trace: Haplotype) -> Option<String> {
        let known_id = self.individuals
            .iter()
            .find(|i| i.haplotype_class == HaplotypeClass::Known)?
            .id
            .clone();

        // Allocate a unique numeric ID.
        let existing: std::collections::HashSet<u64> = self.individuals
            .iter()
            .filter_map(|i| i.id.parse::<u64>().ok())
            .collect();
        let new_id = (0u64..).find(|n| !existing.contains(n)).unwrap_or(0).to_string();
        let name = format!("trace_child_{}", new_id);

        self.add_individual(new_id.clone(), name.clone());
        if let Some(ind) = self.get_individual_by_id_mut(&new_id) {
            ind.haplotype = trace;
            ind.haplotype_class = HaplotypeClass::Suspect;
        }
        self.add_relationship(known_id, new_id);
        Some(name)
    }

    /// Remove an individual and all relationships involving them.
    ///
    /// Does NOT cascade to children — callers are responsible for ordering
    /// (level-order traversal removes parents before children, so orphaned
    /// children will be removed in a subsequent pass if they are also marked).
    pub fn remove_individual(&mut self, id: &str) {
        self.individuals.retain(|i| i.id != id);
        self.relationships.retain(|r| r.parent_id != id && r.child_id != id);
    }

    /// Prune the pedigree by removing irrelevant known individuals.
    ///
    /// Mirrors Python's `remove_irrelevant_individuals()` in models.py:732-793.
    ///
    /// `inside=true`: Removes known/suspect individuals for which EVERY path to
    /// every unknown individual passes through another known individual.
    ///
    /// `inside=false` (extended pedigree): Keeps only nodes on direct paths
    /// from `last_child_name` to any known/suspect individual (paths must pass
    /// only through unknowns).
    ///
    /// Returns the name of the pedigree root to use for subsequent rerooting:
    /// - If the suspect (class == Suspect) survives pruning: returns suspect's name.
    /// - If the suspect is removed: returns the first non-removed known individual's name.
    /// - Returns None if neither exists.
    pub fn remove_irrelevant_individuals(
        &mut self,
        inside: bool,
        last_child_name: Option<&str>,
    ) -> Option<String> {
        use crate::graph::{bfs_layers, shortest_path_undirected};

        let root_ids = self.roots();
        let root_id = root_ids.first()?.to_string();

        let unknown_ids: Vec<String> = self.individuals
            .iter()
            .filter(|i| i.haplotype_class == HaplotypeClass::Unknown)
            .map(|i| i.id.clone())
            .collect();

        let mut to_remove: std::collections::HashSet<String> = std::collections::HashSet::new();

        if inside {
            // For each known/suspect individual: check if it is "irrelevant" —
            // i.e., for every unknown, the shortest undirected path from this
            // known to that unknown passes through at least one other known node.
            let known_suspect_ids: Vec<String> = self.individuals
                .iter()
                .filter(|i| matches!(i.haplotype_class, HaplotypeClass::Known | HaplotypeClass::Suspect))
                .map(|i| i.id.clone())
                .collect();

            // Also remove excluded individuals whose entire subtree is excluded.
            for ind in &self.individuals {
                if ind.exclude {
                    let desc = crate::graph::descendants(self, &ind.id);
                    if desc.iter().all(|did| {
                        self.get_individual_by_id(did)
                            .map(|d| d.exclude)
                            .unwrap_or(true)
                    }) {
                        to_remove.insert(ind.id.clone());
                        for did in desc {
                            to_remove.insert(did);
                        }
                    }
                }
            }

            // Mark irrelevant known individuals.
            // Excluded individuals are never checked for relevance (mirrors Python's elif).
            'outer: for known_id in known_suspect_ids.iter().filter(|id| {
                self.get_individual_by_id(id).map(|i| !i.exclude).unwrap_or(false)
            }) {
                if to_remove.contains(known_id) {
                    continue;
                }
                for unknown_id in &unknown_ids {
                    if let Some(path) = shortest_path_undirected(self, known_id, unknown_id) {
                        // Interior nodes are path[1..path.len()-1]
                        let has_known_interior = path[1..path.len().saturating_sub(1)]
                            .iter()
                            .any(|nid| {
                                self.get_individual_by_id(nid)
                                    .map(|i| i.haplotype_class == HaplotypeClass::Known)
                                    .unwrap_or(false)
                            });
                        if !has_known_interior {
                            // This known individual is directly connected to an unknown —
                            // it is relevant; do not remove.
                            continue 'outer;
                        }
                    }
                }
                // All paths to unknowns are mediated by another known — irrelevant.
                to_remove.insert(known_id.clone());
            }
        } else {
            // Outside mode: keep only nodes on direct-unknown paths from last_child to each known.
            let last_child_id = last_child_name
                .and_then(|name| self.get_individual_by_name(name))
                .map(|i| i.id.clone())?;

            let known_suspect_ids: Vec<String> = self.individuals
                .iter()
                .filter(|i| matches!(i.haplotype_class, HaplotypeClass::Known | HaplotypeClass::Suspect))
                .map(|i| i.id.clone())
                .collect();

            let mut keep: std::collections::HashSet<String> = std::collections::HashSet::new();
            for known_id in &known_suspect_ids {
                if let Some(path) = shortest_path_undirected(self, &last_child_id, known_id) {
                    // Keep only if interior nodes are all unknown.
                    let all_unknown_interior = path[1..path.len().saturating_sub(1)]
                        .iter()
                        .all(|nid| {
                            self.get_individual_by_id(nid)
                                .map(|i| i.haplotype_class == HaplotypeClass::Unknown)
                                .unwrap_or(false)
                        });
                    if all_unknown_interior {
                        for nid in &path {
                            keep.insert(nid.clone());
                        }
                    }
                }
            }

            for ind in &self.individuals {
                if !keep.contains(&ind.id) {
                    to_remove.insert(ind.id.clone());
                }
            }
        }

        // Identify the suspect before removal.
        let suspect_id = self.individuals
            .iter()
            .find(|i| i.haplotype_class == HaplotypeClass::Suspect)
            .map(|i| i.id.clone());
        let suspect_name = suspect_id.as_ref()
            .and_then(|sid| self.get_individual_by_id(sid))
            .map(|i| i.name.clone());

        // If suspect is being removed, find a replacement root.
        let non_removed_known_name = if suspect_id.as_deref().map(|sid| to_remove.contains(sid)).unwrap_or(false) {
            let mut sorted_known: Vec<_> = self.individuals
                .iter()
                .filter(|i| i.haplotype_class == HaplotypeClass::Known && !to_remove.contains(&i.id))
                .map(|i| (i.id.clone(), i.name.clone()))
                .collect();
            sorted_known.sort_by_key(|(id, _)| id.clone());
            sorted_known.into_iter().next().map(|(_, name)| name)
        } else {
            None
        };

        // Remove in BFS (level) order so parents are removed before children.
        let level_order: Vec<String> = bfs_layers(self, &root_id)
            .unwrap_or_default()
            .into_iter()
            .flatten()
            .collect();
        let to_remove_vec: Vec<String> = level_order
            .into_iter()
            .filter(|id| to_remove.contains(id))
            .collect();
        for id in to_remove_vec {
            self.remove_individual(&id);
        }

        // Return root name for subsequent rerooting.
        non_removed_known_name.or(suspect_name)
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
    /// Optional RNG seed for reproducible runs (None = use default deterministic seeds)
    pub seed: Option<u64>,
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
            seed: None,
        }
    }
}

// ---------------------------------------------------------------------------
// StageStats
// ---------------------------------------------------------------------------

/// Per-stage performance statistics collected during simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageStats {
    /// Total iterations run per model (all 3 values are the same: trial_nr * batch_length)
    pub iterations_per_model: Vec<u64>,
    /// Final running mean of each model, formatted as 4-decimal scientific notation
    pub model_probabilities: Vec<String>,
    /// Wall-clock seconds for this stage
    pub runtime_secs: f64,
}

impl StageStats {
    pub fn total_iterations(&self) -> u64 {
        self.iterations_per_model.iter().sum()
    }

    pub fn formatted_runtime(&self) -> String {
        let s = self.runtime_secs;
        if s < 60.0 {
            format!("{:.1}s", s)
        } else {
            format!("{:.0}m {:.0}s", (s / 60.0).floor(), s % 60.0)
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
    /// Per-stage performance stats
    pub pedigree_stats: StageStats,
    pub inside_stats: Option<StageStats>,
    pub extended_pedigree_stats: Option<StageStats>,
    pub outside_stats: Option<StageStats>,
    /// Total wall-clock seconds for the full simulation
    pub total_runtime_secs: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchProbabilities {
    /// match_count → probability
    pub probabilities: HashMap<u32, Decimal>,
    pub average_pedigree_probability: Decimal,
}
