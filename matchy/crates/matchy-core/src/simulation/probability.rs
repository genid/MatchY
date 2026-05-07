/// Edge probability calculation: P(child haplotype | parent haplotype).
///
/// Uses permutation-based matching for multi-copy markers.
/// For single-copy markers this reduces to a simple step probability lookup.
use crate::{Allele, Haplotype, MarkerSet};
use crate::simulation::mutation::neutral_step_probabilities;
use itertools::Itertools; // add to Cargo.toml

/// Probability of observing the child alleles given the parent alleles for a
/// single marker.
///
/// For a marker with n copies, sums over all n! permutations of the parent
/// alleles and for each permutation computes the product of per-copy step
/// probabilities, then takes the average.
///
/// Mirrors `calculate_mutation_probability` in Python's models.py.
pub fn calculate_mutation_probability(
    parent_alleles: &[&Allele],
    child_alleles: &[&Allele],
    mu: f64,
    two_step_fraction: f64,
) -> f64 {
    if parent_alleles.is_empty() || child_alleles.is_empty() {
        return 1.0;
    }
    debug_assert_eq!(
        parent_alleles.len(),
        child_alleles.len(),
        "copy count must match"
    );

    let probs = neutral_step_probabilities(mu, two_step_fraction);
    let n = parent_alleles.len();

    if n == 1 {
        // Fast path: single copy.
        // Intermediate value must match — mirrors Python models.py:1047
        if child_alleles[0].intermediate_value != parent_alleles[0].intermediate_value {
            return 0.0;
        }
        let step = child_alleles[0].value - parent_alleles[0].value;
        return step_probability(step, &probs);
    }

    // Permute CHILD alleles and deduplicate — exact mirror of Python's
    // generate_unique_matchings (models.py:1021-1032).
    // Deduplication matters when two child alleles have the same (value,
    // intermediate_value): Python's `seen` set discards the duplicate permutation,
    // so we must do the same to avoid double-counting.
    let mut seen: std::collections::HashSet<Vec<(i32, i32)>> = std::collections::HashSet::new();
    child_alleles
        .iter()
        .permutations(n)
        .filter(|child_perm| {
            let key: Vec<(i32, i32)> = child_perm
                .iter()
                .map(|a| (a.value, a.intermediate_value.unwrap_or(0)))
                .collect();
            seen.insert(key)
        })
        .map(|child_perm| {
            parent_alleles
                .iter()
                .zip(child_perm.iter())
                .map(|(pa, ca)| {
                    // Intermediate value must match — mirrors Python models.py:1047
                    if pa.intermediate_value != ca.intermediate_value {
                        return 0.0;
                    }
                    let step = ca.value - pa.value;
                    step_probability(step, &probs)
                })
                .product::<f64>()
        })
        .sum()
}

/// Look up the probability for a given mutation step.
/// Steps outside ±2 have probability 0.
fn step_probability(step: i32, probs: &[f64; 5]) -> f64 {
    match step {
        -2 => probs[0],
        -1 => probs[1],
        0 => probs[2],
        1 => probs[3],
        2 => probs[4],
        _ => 0.0,
    }
}

/// Probability of the full child haplotype given the parent haplotype.
///
/// Product of per-marker mutation probabilities.
pub fn get_edge_probability(
    parent: &Haplotype,
    child: &Haplotype,
    marker_set: &MarkerSet,
    two_step_fraction: f64,
) -> f64 {
    let mut probability = 1.0f64;
    for marker in &marker_set.markers {
        let parent_alleles = parent.get_alleles_by_marker_name(&marker.name);
        let child_alleles = child.get_alleles_by_marker_name(&marker.name);

        if parent_alleles.is_empty() || child_alleles.is_empty() {
            continue;
        }
        if parent_alleles.len() != child_alleles.len() {
            // Copy count mismatch — treat as probability 0
            return 0.0;
        }

        let mu = marker.single_copy_mutation_rate();
        let p = calculate_mutation_probability(
            &parent_alleles,
            &child_alleles,
            mu,
            two_step_fraction,
        );
        probability *= p;
    }
    probability
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Allele, Marker};

    fn allele(marker: &Marker, value: i32) -> Allele {
        Allele::new(marker.clone(), value)
    }

    /// Failing case: only valid permutation requires a ±2 step.
    /// With tsf=0.03 this should be non-zero; with tsf=0 it should be zero.
    #[test]
    fn multi_copy_requires_two_step() {
        let m = Marker::new("TEST", 0.002).with_copies(3);
        let mu = m.single_copy_mutation_rate();

        let pa: Vec<Allele> = [59, 63, 68].iter().map(|&v| allele(&m, v)).collect();
        let ca: Vec<Allele> = [60, 61, 67].iter().map(|&v| allele(&m, v)).collect();
        let par: Vec<&Allele> = pa.iter().collect();
        let chi: Vec<&Allele> = ca.iter().collect();

        // With tsf=0.03: only valid permutation is (60,61,67) → steps (+1,-2,-1)
        // p(-2) = mu*0.03/2 > 0, so result should be non-zero.
        let p_with_tsf = calculate_mutation_probability(&par, &chi, mu, 0.03);
        assert!(p_with_tsf > 0.0, "expected non-zero with tsf=0.03, got {}", p_with_tsf);

        // With tsf=0: p(±2)=0, so result should be zero.
        let p_no_tsf = calculate_mutation_probability(&par, &chi, mu, 0.0);
        assert_eq!(p_no_tsf, 0.0, "expected zero with tsf=0.0, got {}", p_no_tsf);
    }

    #[test]
    fn multi_copy_no_two_step_needed() {
        let m = Marker::new("TEST", 0.002).with_copies(2);
        let mu = m.single_copy_mutation_rate();

        // parent=[37,39] child=[40,40]: best steps are (+3,+1) → always zero
        let pa: Vec<Allele> = [37, 39].iter().map(|&v| allele(&m, v)).collect();
        let ca: Vec<Allele> = [40, 40].iter().map(|&v| allele(&m, v)).collect();
        let par: Vec<&Allele> = pa.iter().collect();
        let chi: Vec<&Allele> = ca.iter().collect();

        let p = calculate_mutation_probability(&par, &chi, mu, 0.03);
        assert_eq!(p, 0.0, "expected zero (>2 step required), got {}", p);
    }

    #[test]
    fn multi_copy_all_within_one_step() {
        let m = Marker::new("TEST", 0.002).with_copies(2);
        let mu = m.single_copy_mutation_rate();

        // parent=[37,39] child=[37,40]: best perm (37,40): steps (0,+1) → non-zero
        let pa: Vec<Allele> = [37, 39].iter().map(|&v| allele(&m, v)).collect();
        let ca: Vec<Allele> = [37, 40].iter().map(|&v| allele(&m, v)).collect();
        let par: Vec<&Allele> = pa.iter().collect();
        let chi: Vec<&Allele> = ca.iter().collect();

        let p = calculate_mutation_probability(&par, &chi, mu, 0.0);
        assert!(p > 0.0, "expected non-zero (only ±1 steps needed), got {}", p);
    }
}

/// Compute edge probabilities for all parent→child relationships in the pedigree.
/// Returns a map: (parent_id, child_id) → probability.
pub fn get_all_edge_probabilities(
    pedigree: &crate::Pedigree,
    marker_set: &MarkerSet,
    two_step_fraction: f64,
) -> std::collections::HashMap<(String, String), f64> {
    pedigree
        .relationships
        .iter()
        .filter_map(|rel| {
            let parent = pedigree.get_individual_by_id(&rel.parent_id)?;
            let child = pedigree.get_individual_by_id(&rel.child_id)?;
            let p = get_edge_probability(
                &parent.haplotype,
                &child.haplotype,
                marker_set,
                two_step_fraction,
            );
            Some(((rel.parent_id.clone(), rel.child_id.clone()), p))
        })
        .collect()
}
