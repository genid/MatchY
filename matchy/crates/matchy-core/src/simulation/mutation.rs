/// Mutation probability tables and haplotype mutation sampling.
///
/// Mutation model: step changes in allele values.
/// Steps = {-2, -1, 0, +1, +2}
/// p(0)  = 1 - mu
/// p(±1) = (mu * (1 - two_step_fraction)) / 2
/// p(±2) = (mu * two_step_fraction) / 2
///
/// Bias model (matching Python): adds target_mass to the chosen direction,
/// subtracts from the no-step (0) probability. The distribution remains
/// unnormalized; WeightedIndex handles that.
use crate::{Allele, Bias, BiasDirection, Haplotype};
use rand::distributions::WeightedIndex;
use rand::prelude::*;
use std::collections::HashMap;

/// The five possible mutation steps.
pub const STEPS: [i32; 5] = [-2, -1, 0, 1, 2];

/// Compute neutral (unbiased) step probabilities for a single-copy marker.
///
/// Returns weights for STEPS = [-2, -1, 0, +1, +2].
pub fn neutral_step_probabilities(mu: f64, two_step_fraction: f64) -> [f64; 5] {
    let p_no_mutation = (1.0 - mu).max(0.0);
    let p_one_step = (mu * (1.0 - two_step_fraction) / 2.0).max(0.0);
    let p_two_step = (mu * two_step_fraction / 2.0).max(0.0);
    [p_two_step, p_one_step, p_no_mutation, p_one_step, p_two_step]
}

/// Apply a Bias to the step probability distribution.
///
/// Matches Python's approach: add target_mass to the biased direction's
/// one-step probability, subtract from the zero-step probability.
/// The resulting distribution is unnormalized (valid for WeightedIndex).
pub fn apply_bias(probs: &[f64; 5], bias: &Bias) -> [f64; 5] {
    let mut adjusted = *probs;
    let target = bias.target_mass;
    match bias.direction {
        BiasDirection::Up => {
            // Shift target_mass from no-step to +1 step
            adjusted[3] += target; // p(+1)
            adjusted[2] = (adjusted[2] - target).max(0.0); // p(0)
        }
        BiasDirection::Down => {
            // Shift target_mass from no-step to -1 step
            adjusted[1] += target; // p(-1)
            adjusted[2] = (adjusted[2] - target).max(0.0); // p(0)
        }
    }
    adjusted
}

/// Mutate a single allele by sampling a step from the distribution.
///
/// Returns (new_allele, weighted_prob, unweighted_prob).
/// - weighted_prob: probability under the biased distribution (for importance sampling)
/// - unweighted_prob: probability under the neutral distribution
///
/// Intermediate values are preserved only for zero-step mutations.
pub fn mutate_allele<R: Rng>(
    allele: &Allele,
    copy_nr: usize,
    mu: f64,
    two_step_fraction: f64,
    biases: &[Bias],
    rng: &mut R,
) -> (Allele, f64, f64) {
    let neutral = neutral_step_probabilities(mu, two_step_fraction);

    // Apply biases for this marker and copy number (matching Python)
    let biased = biases
        .iter()
        .filter(|b| b.marker.name == allele.marker.name && b.copy_nr as usize == copy_nr)
        .fold(neutral, |acc, b| apply_bias(&acc, b));

    let dist = WeightedIndex::new(&biased).expect("valid step probabilities");
    let step_idx = dist.sample(rng);
    let step = STEPS[step_idx];

    let weighted_prob = biased[step_idx];
    let unweighted_prob = neutral[step_idx];

    let new_value = (allele.value + step).max(1); // allele values are ≥ 1
    // Always preserve intermediate_value regardless of step — mirrors Python
    // simulation.py:68: Allele(marker, mutated_value, source_allele.intermediate_value)
    let new_intermediate = allele.intermediate_value;

    let new_allele = Allele {
        marker: allele.marker.clone(),
        value: new_value,
        intermediate_value: new_intermediate,
        mutation_value: Some(step),
        mutation_probability: Some(weighted_prob),
    };

    (new_allele, weighted_prob, unweighted_prob)
}

/// Mutate all alleles in a haplotype by one generation using importance sampling.
///
/// For each marker, mutates each copy independently with bias.
///
/// Returns (new_haplotype, w_edge, u_edge) where:
/// - w_edge = product of biased step probabilities
/// - u_edge = product of neutral step probabilities
///
/// The importance weight factor = u_edge / w_edge (reciprocal of IS correction).
pub fn mutate_haplotype<R: Rng>(
    haplotype: &Haplotype,
    marker_set: &crate::MarkerSet,
    two_step_fraction: f64,
    biases: &[Bias],
    rng: &mut R,
) -> (Haplotype, f64, f64) {
    let mut new_haplotype = Haplotype::new();
    let mut w_edge = 1.0f64; // weighted (biased) probability product
    let mut u_edge = 1.0f64; // unweighted (neutral) probability product

    for marker in &marker_set.markers {
        let marker_name = &marker.name;
        // Use sorted allele order to match bias copy_nr assignment in get_biases_for_individual.
        let alleles = haplotype.get_alleles_by_marker_name(marker_name);
        if alleles.is_empty() {
            continue;
        }

        let mu = marker.single_copy_mutation_rate();

        for (copy_nr, allele) in alleles.iter().enumerate() {
            let (new_allele, w, u) = mutate_allele(
                allele,
                copy_nr,
                mu,
                two_step_fraction,
                biases,
                rng,
            );
            w_edge *= w;
            u_edge *= u;
            new_haplotype
                .alleles
                .entry(marker_name.clone())
                .or_default()
                .push(new_allele);
        }
    }

    (new_haplotype, w_edge, u_edge)
}

/// Mutate all alleles using pre-computed neutral step probability tables.
///
/// Identical to `mutate_haplotype` but skips recomputing `neutral_step_probabilities`
/// on every call. Use inside a simulation batch where the tables are precomputed
/// once and reused across all iterations.
pub fn mutate_haplotype_precomputed<R: Rng>(
    haplotype: &Haplotype,
    marker_set: &crate::MarkerSet,
    neutral_tables: &HashMap<String, [f64; 5]>,
    biases: &[Bias],
    rng: &mut R,
) -> (Haplotype, f64, f64) {
    let mut new_haplotype = Haplotype::new();
    let mut w_edge = 1.0f64;
    let mut u_edge = 1.0f64;

    for marker in &marker_set.markers {
        let marker_name = &marker.name;
        // Use sorted allele order to match bias copy_nr assignment in get_biases_for_individual.
        let alleles = haplotype.get_alleles_by_marker_name(marker_name);
        if alleles.is_empty() {
            continue;
        }
        let neutral = match neutral_tables.get(marker_name) {
            Some(n) => n,
            None => continue,
        };

        for (copy_nr, allele) in alleles.iter().enumerate() {
            let biased = biases
                .iter()
                .filter(|b| b.marker.name == *marker_name && b.copy_nr as usize == copy_nr)
                .fold(*neutral, |acc, b| apply_bias(&acc, b));

            let dist = WeightedIndex::new(&biased).expect("valid step probabilities");
            let step_idx = dist.sample(rng);
            let step = STEPS[step_idx];
            let weighted_prob = biased[step_idx];
            let unweighted_prob = neutral[step_idx];

            let new_allele = Allele {
                marker: allele.marker.clone(),
                value: (allele.value + step).max(1),
                intermediate_value: allele.intermediate_value,
                mutation_value: Some(step),
                mutation_probability: Some(weighted_prob),
            };

            w_edge *= weighted_prob;
            u_edge *= unweighted_prob;
            new_haplotype
                .alleles
                .entry(marker_name.clone())
                .or_default()
                .push(new_allele);
        }
    }

    (new_haplotype, w_edge, u_edge)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Marker, Haplotype};
    use rand::SeedableRng;
    use rand::rngs::SmallRng;

    /// Regression test: bias copy_nr must align with sorted allele order.
    ///
    /// Constructs a 2-copy marker haplotype with alleles stored in REVERSE sorted
    /// order (high value first, low value second — simulating value drift across
    /// generations that reorders copies). A bias for copy_nr=0 ("Up") targets the
    /// sorted-first (lowest-value) allele. Before this fix, the bias was applied to
    /// the first-stored allele (the high-value one) instead.
    ///
    /// Verifies: the biased copy (sorted copy 0, lowest value) drifts upward on
    /// average while the unbiased copy (sorted copy 1) stays flat.
    #[test]
    fn bias_applied_to_sorted_copy_not_storage_copy() {
        let m = Marker::new("TEST", 0.002).with_copies(2);
        let tsf = 0.03;

        // Build haplotype with alleles inserted in reverse sorted order:
        // storage: [40, 30]  but sorted order: [30, 40]
        let mut hap = Haplotype::new();
        hap.add_allele(m.clone(), 40, None); // storage copy 0 (high value)
        hap.add_allele(m.clone(), 30, None); // storage copy 1 (low value)

        let marker_set = crate::MarkerSet { markers: vec![m.clone()] };
        let neutral_tables: HashMap<String, [f64; 5]> = [(
            "TEST".to_string(),
            neutral_step_probabilities(m.single_copy_mutation_rate(), tsf),
        )].into_iter().collect();

        // Bias copy_nr=0 Up — should target sorted copy 0 = value 30.
        let biases = vec![Bias {
            marker: m.clone(),
            copy_nr: 0,
            direction: BiasDirection::Up,
            target_mass: 0.4,
        }];

        let n = 2000u32;
        let mut sum_low = 0i64;
        let mut sum_high = 0i64;
        let mut rng = SmallRng::seed_from_u64(42);

        for _ in 0..n {
            let (new_hap, _, _) = mutate_haplotype_precomputed(
                &hap, &marker_set, &neutral_tables, &biases, &mut rng,
            );
            let alleles = new_hap.get_alleles_by_marker_name("TEST");
            sum_low  += alleles[0].value as i64 - 30;
            sum_high += alleles[1].value as i64 - 40;
        }

        let avg_low  = sum_low  as f64 / n as f64;
        let avg_high = sum_high as f64 / n as f64;

        // Biased copy (sorted 0, value 30) must drift positive.
        assert!(avg_low > 0.05,
            "sorted copy 0 (low, biased Up) should drift positive; got {:.4}", avg_low);
        // Unbiased copy (sorted 1, value 40) must stay near zero.
        assert!(avg_high.abs() < 0.05,
            "sorted copy 1 (high, unbiased) should stay flat; got {:.4}", avg_high);
    }
}
