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
use crate::{Allele, Bias, BiasDirection, Haplotype, Marker};
use rand::distributions::WeightedIndex;
use rand::prelude::*;

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
        let alleles = match haplotype.alleles.get(marker_name) {
            Some(a) => a,
            None => continue,
        };

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
