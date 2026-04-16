/// Hungarian algorithm wrapper for optimal allele matching.
///
/// Mirrors Python's SimpleDifferenceMatrix which uses
/// scipy.optimize.linear_sum_assignment to find the optimal pairing of
/// parent alleles to child alleles, then computes the signed mutation
/// direction for each pair.
///
/// INTERMEDIATE_MISMATCH_PENALTY inflates the cost for pairs where the
/// intermediate values differ, steering the assignment away from them.
/// After the assignment, the actual mutation direction is computed from
/// the raw allele values (ignoring the penalty).
use crate::{Allele, Bias, BiasDirection, Marker};
use pathfinding::prelude::{kuhn_munkres_min, Matrix};

const INTERMEDIATE_MISMATCH_PENALTY: i32 = 1000;

/// Represents the difference matrix between parent and child alleles.
pub struct SimpleDifferenceMatrix<'a> {
    pub parent_alleles: Vec<&'a Allele>,
    pub child_alleles: Vec<&'a Allele>,
}

impl<'a> SimpleDifferenceMatrix<'a> {
    pub fn new(parent_alleles: Vec<&'a Allele>, child_alleles: Vec<&'a Allele>) -> Self {
        Self {
            parent_alleles,
            child_alleles,
        }
    }

    /// Build cost matrix: cost[i][j] = abs(parent[i].value - child[j].value)
    /// + INTERMEDIATE_MISMATCH_PENALTY if intermediate values differ.
    fn cost_matrix(&self) -> Vec<Vec<i32>> {
        self.parent_alleles
            .iter()
            .map(|p| {
                self.child_alleles
                    .iter()
                    .map(|c| {
                        let base_cost = (p.value - c.value).abs();
                        let penalty = if p.intermediate_value != c.intermediate_value {
                            INTERMEDIATE_MISMATCH_PENALTY
                        } else {
                            0
                        };
                        base_cost + penalty
                    })
                    .collect()
            })
            .collect()
    }

    /// Compute the optimal assignment using the Hungarian algorithm and
    /// return a list of signed mutation values (child.value - parent.value)
    /// plus an intermediate mismatch flag for each pair.
    pub fn calculate_mutations(&self) -> Vec<MutationDelta> {
        if self.parent_alleles.is_empty() || self.child_alleles.is_empty() {
            return Vec::new();
        }

        let cost_rows = self.cost_matrix();
        let nrows = cost_rows.len();
        let ncols = cost_rows[0].len();
        let flat: Vec<i32> = cost_rows.into_iter().flatten().collect();
        let matrix = Matrix::from_vec(nrows, ncols, flat)
            .expect("valid rectangular cost matrix");
        let (_, assignment) = kuhn_munkres_min(&matrix);

        assignment
            .iter()
            .enumerate()
            .map(|(parent_idx, &child_idx)| {
                let parent = self.parent_alleles[parent_idx];
                let child = self.child_alleles[child_idx];
                let step = child.value - parent.value;
                let intermediate_mismatch = parent.intermediate_value != child.intermediate_value;
                MutationDelta {
                    step,
                    intermediate_mismatch,
                    parent_copy: parent_idx as u32,
                    child_copy: child_idx as u32,
                }
            })
            .collect()
    }

    /// Derive Bias structs from the optimal assignment for use in importance sampling.
    /// Returns one Bias per copy of this marker.
    pub fn get_biases(&self, marker: &Marker) -> Vec<Bias> {
        let deltas = self.calculate_mutations();
        deltas
            .iter()
            .enumerate()
            .filter_map(|(copy_nr, delta)| {
                if delta.step == 0 {
                    return None; // no bias needed for no-mutation
                }
                let direction = if delta.step > 0 {
                    BiasDirection::Up
                } else {
                    BiasDirection::Down
                };
                Some(Bias {
                    marker: marker.clone(),
                    copy_nr: copy_nr as u32,
                    direction,
                    target_mass: 0.0, // filled in by bias.rs
                })
            })
            .collect()
    }
}

/// The signed number of steps for one parent→child allele pair.
#[derive(Debug, Clone)]
pub struct MutationDelta {
    pub step: i32,
    pub intermediate_mismatch: bool,
    pub parent_copy: u32,
    pub child_copy: u32,
}
