# MatchY — Implementation Reference

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Repository Layout](#2-repository-layout)
3. [Core Data Structures (`matchy-core/models`)](#3-core-data-structures)
4. [Pedigree Graph (`matchy-core/graph`)](#4-pedigree-graph)
5. [Hungarian Assignment (`matchy-core/hungarian`)](#5-hungarian-assignment)
6. [Mutation Model (`simulation/mutation`)](#6-mutation-model)
7. [Edge Probability (`simulation/probability`)](#7-edge-probability)
8. [Importance Sampling Bias (`simulation/bias`)](#8-importance-sampling-bias)
9. [Monte Carlo Engine (`simulation/importance`)](#9-monte-carlo-engine)
10. [Three-Model Ensemble & Convergence (`simulation/convergence`)](#10-three-model-ensemble--convergence)
11. [Top-Level Simulation Flow (`simulation/mod`)](#11-top-level-simulation-flow)
12. [File I/O (`matchy-io`)](#12-file-io)
13. [Report Generation (`matchy-report`)](#13-report-generation)
14. [CLI (`matchy-cli`)](#14-cli)
15. [Tauri Desktop App — Rust Backend (`matchy-app`)](#15-tauri-desktop-app--rust-backend)
16. [Frontend — React / TypeScript (`matchy-app/frontend`)](#16-frontend--react--typescript)

---

## 1. Project Overview

MatchY is a forensic/genealogical Y-STR kinship probability tool. Given a pedigree of related individuals (some with known Y-STR haplotypes, some unknown), it estimates the probability that an unknown donor within the pedigree matches a given profile — the **suspect** or **trace** haplotype.

The key scientific question answered: *What is the probability that a randomly-selected unknown male in the pedigree produced the observed Y-STR profile?*

Two main outputs are computed:

- **Inside match probability** — the probability that the suspect haplotype matches one or more unknowns in the pedigree, broken down by number of matches and per individual.
- **Outside match probability** — an analogous probability calculated for an "extended" pedigree that includes one extra generation, representing the possibility that the true donor is not captured by the known pedigree.

The engine uses **Monte Carlo importance sampling** with a **three-model ensemble** and adaptive convergence checking. Computations run in parallel via Rayon.

The project ships as:
- A Rust library crate (`matchy-core`) containing the entire simulation engine.
- A CLI tool (`matchy-cli`) driven by a TOML/INI config file.
- A Tauri desktop application (`matchy-app`) with a React/TypeScript frontend.

---

## 2. Repository Layout

```
matchy/
├── Cargo.toml                      # workspace (edition 2021, MIT)
└── crates/
    ├── matchy-core/                # simulation engine, data models, graph
    │   └── src/
    │       ├── lib.rs
    │       ├── models.rs
    │       ├── graph.rs
    │       ├── hungarian.rs
    │       └── simulation/
    │           ├── mod.rs          # top-level entry point
    │           ├── bias.rs         # importance-sampling bias
    │           ├── convergence.rs  # three-model ensemble & convergence
    │           ├── importance.rs   # per-batch Monte Carlo loops
    │           ├── mutation.rs     # mutation probability model
    │           └── probability.rs  # edge probability calculation
    ├── matchy-io/                  # file format parsers, embedded kits/rates
    ├── matchy-report/              # Minijinja HTML report rendering
    ├── matchy-cli/                 # CLI entry point
    └── matchy-app/
        ├── src/                    # Tauri Rust backend
        └── frontend/               # React 18 + TypeScript frontend
```

Key dependencies: `rayon` (parallelism), `rand`/`rand_distr` (sampling), `petgraph` (DAG), `pathfinding` (Hungarian algorithm), `rust_decimal` (exact decimal arithmetic), `serde`/`serde_json` (serialisation), `minijinja` (report templates), `lru`/`rustc-hash` (caching).

---

## 3. Core Data Structures

### 3.1 `Marker`

```rust
pub struct Marker {
    pub name: String,
    pub mutation_rate: f64,          // per-generation, all-copy rate
    pub number_of_copies: Option<u32>,
}
```

When a marker has multiple copies (e.g. `DYF387S1` with 2 copies), the all-copy mutation rate must be decomposed into a per-copy rate:

```rust
pub fn single_copy_mutation_rate(&self) -> f64 {
    match self.number_of_copies {
        Some(n) if n > 1 => 1.0 - (1.0 - self.mutation_rate).powf(1.0 / n as f64),
        _ => self.mutation_rate,
    }
}
```

**Pseudocode:**
```
per_copy_rate = 1 - (1 - all_copy_rate)^(1 / n_copies)
```

The relationship is derived from: the probability of *no* mutation in *any* of `n` independent copies is `(1 - per_copy_rate)^n`, and this must equal `(1 - all_copy_rate)`.

### 3.2 `Allele`

```rust
pub struct Allele {
    pub marker: Marker,
    pub value: i32,                       // integer repeat count
    pub intermediate_value: Option<i32>,  // e.g. "15.2" → value=15, intermediate=2
    pub mutation_value: Option<i32>,      // set during simulation
    pub mutation_probability: Option<f64>,
}
```

Alleles are compared and ordered by `(value, intermediate_value_or_zero)`. Display format: `"15"` for simple, `"15.2"` for intermediate alleles.

### 3.3 `Haplotype`

```rust
pub struct Haplotype {
    pub alleles: HashMap<String, Vec<Allele>>,  // marker_name → sorted Vec<Allele>
}
```

Key methods:

**`allelic_difference(other) -> Option<i32>`**  
The sum of absolute differences between paired alleles across all markers (after sorting alleles within each marker). Returns `None` if marker sets are incompatible.

```
distance = Σ_markers Σ_copies |self.allele[i].value - other.allele[i].value|
```

**Hashing** is deterministic: markers are sorted lexicographically, alleles within each marker are sorted by `(value, intermediate_norm)`, so the same logical haplotype always produces the same hash regardless of insertion order.

### 3.4 `Individual`

```rust
pub struct Individual {
    pub id: String,
    pub name: String,
    pub haplotype: Haplotype,
    pub haplotype_class: HaplotypeClass,  // Unknown | Known | Suspect | Estimated | Fixed | Excluded
    pub exclude: bool,
    pub picking_probability: Option<Decimal>,
    pub closest_known_individuals: Vec<String>,
    pub closest_known_distance: Option<u32>,
}
```

An individual is **known** (has an observed haplotype) if `haplotype_class` is `Known`, `Suspect`, or `Fixed`. This distinction matters throughout the simulation — only unknowns are simulated; known individuals contribute fixed probabilities via edge calculations.

### 3.5 `Pedigree`

```rust
pub struct Pedigree {
    pub individuals: Vec<Individual>,
    pub relationships: Vec<Relationship>,     // directed parent → child edges
    pub picking_probabilities: HashMap<String, Decimal>,
    pub trace_haplotype: Option<Haplotype>,
}
```

#### `extend_pedigree() -> String`

Extends the pedigree upward and outward to model the possibility that the donor is outside the known pedigree. Returns the name of the new "outside" leaf node.

```
Algorithm:
1. Find current root(s) via roots().
2. BFS from root; find the shallowest layer containing a Known individual → first_known_depth.
3. target_depth = first_known_depth + 1
4. Add a new root ("new_root_N") above the current root.
5. Add a chain of (target_depth - 1) unknown children below the new root,
   down to the same level as the original root.
6. Return the name of the last added child ("outside" individual).
```

#### `reroot(new_root_id)`

Re-orients all edges to flow away from `new_root_id`. This is needed after pruning so that BFS traversal during simulation starts from the correct root.

```
Algorithm:
1. DFS from new_root_id treating all relationships as undirected.
2. Reconstruct relationships: for each traversal edge (u → v), add relationship parent=u, child=v.
3. Demote prior Suspect individuals back to Known.
4. Promote new_root_id to Suspect (if it was Known).
```

#### `calculate_picking_probabilities(trace_mode)`

Determines how likely each unknown individual is to be the true donor, used as a prior in the Monte Carlo sampling.

**Trace mode** (uniform prior):
```
P(individual) = 1 / n_unknowns
```

**Suspect mode** (similarity-weighted prior):
```
Algorithm (mirrors original Python models.py):
1. For each unknown U: find shortest undirected path from suspect to U.
2. Copy haplotypes along that path (Estimated class) to give U a rough
   expected haplotype based on proximity to the suspect.
3. For each unknown U:
   a. Fix U's haplotype to the suspect haplotype (class = Fixed).
   b. Compute sum of allelic distances of every other unknown to its
      path-parent. (Add 1 to avoid zeros.)
   c. distance_score(U) = sum_of_distances_when_U_is_fixed + 1
4. Normalise:
   raw_prob(U) = (max_score - score(U) + 1) / (max_score + 1)
   P(U) = raw_prob(U) / Σ raw_prob
```

The intuition: unknowns that are genetically "close" to the suspect (in terms of allelic distance along the pedigree path) get higher picking probabilities.

#### `remove_irrelevant_individuals(inside, last_child_name)`

Prunes the pedigree to remove individuals that contribute no information to the simulation.

**Inside pruning** (`inside=true`):
```
For each Known/Suspect individual K:
    If ALL paths from K to all unknowns pass through at least one other Known individual:
        K is irrelevant → remove K (it is "shielded" from unknowns by other Knowns).
Also remove excluded individuals whose entire subtree is excluded.
Removal order: BFS level order (remove parents before children to avoid orphans).
```

**Outside pruning** (`inside=false`):
```
Keep only nodes that lie on a direct path (all-unknown-interior) from
last_child_name to any Known/Suspect individual.
```

---

## 4. Pedigree Graph

`matchy-core/graph.rs` wraps `petgraph::DiGraph<String, ()>` for DAG operations.

| Function | Description |
|---|---|
| `validate_dag(ped)` | Checks for cycles using petgraph `is_cyclic_directed`. |
| `topological_order(ped)` | Returns individuals in ancestor-before-descendant order. |
| `bfs_layers(ped, root_id)` | Layer 0 = root, layer k = individuals k steps below root. |
| `shortest_path_length(ped, from, to)` | A* on the directed graph. |
| `ancestors(ped, id)` | All ancestors via DFS on reversed edges (inclusive). |
| `descendants(ped, id)` | All descendants via DFS on forward edges (inclusive). |
| `most_recent_common_ancestor(ped, a, b)` | Intersection of ancestor sets; selects deepest. |
| `shortest_path_undirected(ped, from, to)` | BFS on undirected adjacency; returns full path. |
| `is_connected(ped)` | Builds `UnGraph`, checks `connected_components <= 1`. |

---

## 5. Hungarian Assignment

`matchy-core/hungarian.rs` finds the optimal allele pairing between a parent and child for multi-copy markers.

```rust
pub struct SimpleDifferenceMatrix<'a> {
    pub parent_alleles: Vec<&'a Allele>,
    pub child_alleles: Vec<&'a Allele>,
}
```

**Cost matrix:**
```
cost[i][j] = |parent[i].value - child[j].value|
             + INTERMEDIATE_MISMATCH_PENALTY   (= 1000)
               if parent[i].intermediate_value ≠ child[j].intermediate_value
```

The large intermediate mismatch penalty effectively forces the Hungarian algorithm to prefer pairings that preserve intermediate allele consistency, treating intermediate mismatches as highly undesirable.

**`calculate_mutations() -> Vec<MutationDelta>`**  
Runs Kuhn-Munkres (via `pathfinding::kuhn_munkres_min`) on the cost matrix. Returns the optimal assignment as signed mutation steps:
```rust
pub struct MutationDelta {
    pub step: i32,              // child.value - parent.value
    pub intermediate_mismatch: bool,
    pub parent_copy: u32,
    pub child_copy: u32,
}
```

**`get_biases(marker) -> Vec<Bias>`**  
Converts the optimal assignment into directional bias hints (`Up` / `Down`) for the importance sampler.

---

## 6. Mutation Model

`simulation/mutation.rs` implements a **stepwise mutation model** for Y-STR loci. Mutations happen in integer steps of ±1 or ±2 repeat units per generation.

### 6.1 Neutral Step Probabilities

```rust
fn neutral_step_probabilities(mu: f64, two_step_fraction: f64) -> [f64; 5]
// index layout: [-2, -1, 0, +1, +2]
```

```
p(step =  0)  =  1 - mu
p(step = ±1)  =  mu × (1 - two_step_fraction) / 2
p(step = ±2)  =  mu × two_step_fraction / 2
```

**Pseudocode:**
```
neutral_probs[-2] = mu * two_step_fraction / 2
neutral_probs[-1] = mu * (1 - two_step_fraction) / 2
neutral_probs[ 0] = 1 - mu
neutral_probs[+1] = mu * (1 - two_step_fraction) / 2
neutral_probs[+2] = mu * two_step_fraction / 2
```

### 6.2 Applying Importance Bias

To guide the simulation toward configurations consistent with observed data (importance sampling), the neutral distribution is modified:

```rust
fn apply_bias(probs: &[f64; 5], bias: &Bias) -> [f64; 5]
```

```
if bias.direction == Up:
    probs[+1] += bias.target_mass
    probs[ 0]  = max(0, probs[0] - bias.target_mass)
if bias.direction == Down:
    probs[-1] += bias.target_mass
    probs[ 0]  = max(0, probs[0] - bias.target_mass)
```

The result is an unnormalised distribution; `rand::distributions::WeightedIndex` handles the normalisation implicitly during sampling.

### 6.3 Mutating a Single Allele

```rust
fn mutate_allele(allele, copy_nr, mu, two_step_fraction, biases, rng)
    -> (Allele, f64, f64)   // (new_allele, biased_prob, neutral_prob)
```

```
Algorithm:
1. Compute neutral_probs(mu, two_step_fraction).
2. For each Bias b that applies to this copy_nr:
       probs = apply_bias(probs, b)
3. Sample step index ~ WeightedIndex(probs).
4. new_value = max(1, allele.value + step)       // floor at 1
5. weighted_prob  = probs[step_idx] / sum(probs)
6. neutral_prob   = neutral_probs[step_idx] / sum(neutral_probs)
7. Return (new_allele with new_value, weighted_prob, neutral_prob)
```

Intermediate values are **preserved unchanged** regardless of the step taken. This matches the original Python behaviour: the integer repeat count mutates stepwise, but the fractional micro-variant is inherited as-is.

### 6.4 Mutating a Full Haplotype

```rust
fn mutate_haplotype(haplotype, marker_set, two_step_fraction, biases, rng)
    -> (Haplotype, f64, f64)   // (new_haplotype, w_edge, u_edge)
```

```
w_edge = 1.0   (product of all biased mutation probs)
u_edge = 1.0   (product of all neutral mutation probs)

For each marker in marker_set:
    mu_copy = marker.single_copy_mutation_rate()
    For each copy i of this marker:
        (new_allele, w, u) = mutate_allele(allele[i], i, mu_copy, ...)
        w_edge *= w
        u_edge *= u

Return (new_haplotype, w_edge, u_edge)
```

The ratio `u_edge / w_edge` is the **importance weight** for this edge — how much less (or more) likely this transition was under the neutral model than under the biased model. It corrects for the bias when accumulating statistics.

---

## 7. Edge Probability

`simulation/probability.rs` computes the exact probability of observing a child's haplotype given a parent's haplotype and a marker set.

### 7.1 Single-Copy Marker

```
p(parent → child) = step_probability(child.value - parent.value, neutral_probs)
                    × 0   if intermediate values differ
```

Where `step_probability` looks up the neutral probability for steps in `{-2,-1,0,+1,+2}` and returns 0 for any other step size.

### 7.2 Multi-Copy Marker

For a marker with `n` copies, the child and parent each have `n` alleles. The probability sums over all valid pairings (permutations):

```rust
fn calculate_mutation_probability(parent_alleles, child_alleles, mu, two_step_fraction) -> f64
```

```
Algorithm:
1. Enumerate all n! permutations of child_alleles.
2. Deduplicate permutations (use a HashSet to avoid counting identical
   permutations from repeated allele values).
3. For each unique permutation σ:
       p_perm = Π_i step_probability(child[σ(i)].value - parent[i].value)
                × 0  if any intermediate values differ in the pairing
       total_prob += p_perm
4. Return total_prob
```

**Pseudocode:**
```
for each unique permutation σ of child_alleles:
    p = 1.0
    for i in 0..n:
        step = child_allele[σ(i)].value - parent_allele[i].value
        p *= neutral_step_prob(step, mu, two_step_fraction)
        if child_allele[σ(i)].intermediate ≠ parent_allele[i].intermediate:
            p = 0; break
    total_prob += p
return total_prob
```

### 7.3 Full Haplotype Edge Probability

```
get_edge_probability(parent, child, marker_set, two_step_fraction)
    = Π_markers  calculate_mutation_probability(parent_alleles[m], child_alleles[m], ...)
    = 0  if copy counts mismatch for any marker
```

---

## 8. Importance Sampling Bias

`simulation/bias.rs` provides the adaptive importance-sampling bias infrastructure.

### 8.1 Default Bias Target Mass

```rust
pub fn default_bias_target_mass(distance_to_mrca: u32) -> f64 {
    f64::min(f64::max(0.1, 0.8 / (1.0 + distance_to_mrca as f64)), 0.4)
}
```

**Pseudocode:**
```
target_mass = clamp(0.8 / (1 + distance_to_mrca), min=0.1, max=0.4)
```

Intuition: the further away the most recent common ancestor, the more steps of random drift can blur the signal, so a smaller bias is used to avoid over-constraining the sampler.

### 8.2 Computing Biases for an Individual

`get_biases_for_individual(pedigree, individual_id, marker_set, current_haplotypes, bias_value, fixed_id) -> Vec<Bias>`

```
Algorithm:
1. If bias_value <= 0: return [] (bias disabled).
2. Find this individual's parent in the pedigree.
3. Find all "known" descendants of this individual
   (class Known, Suspect, Fixed, or == fixed_id).
4. Compute MRCA of known descendants; compute distance_to_mrca.
5. For each marker m:
   For each copy c:
     Collect the optimal mutation direction (Up/Down) from parent→known_descendant
     using SimpleDifferenceMatrix.
     If ALL known descendants agree on direction:
         Create Bias(marker=m, copy_nr=c, direction=agreed_dir,
                     target_mass = bias_value ?? default_bias_target_mass(distance_to_mrca))
6. Return list of Bias objects.
```

The rationale: if all downstream known individuals have allele values **above** the parent, the true (unknown) path almost certainly went "Up" at some point. The sampler is therefore biased toward upward steps, then corrected by the importance weight.

### 8.3 Adaptive Bias Schedule

```rust
pub struct AdaptiveBiasSchedule {
    pub trial: u32,
    pub model_biases: [f64; 3],   // start at [0.10, 0.10, 0.10]
}
```

After each convergence trial, the three models are ranked by their last mean probability:
```
Algorithm:
1. Rank models by last_mean (descending).
2. Best model  → bias = 0.05  (lightest bias, trust the sampler more)
3. Middle model → bias = 0.15
4. Worst model  → bias = 0.25  (strongest bias, force more exploration)
```

This ensures the ensemble maintains model diversity: models that are converging well reduce their bias, while poorly-converging models increase it to explore different configurations.

---

## 9. Monte Carlo Engine

`simulation/importance.rs` is the innermost computational loop.

### 9.1 `BatchResult`

Accumulates running statistics for one batch:

```rust
pub struct BatchResult {
    pub weighted_sum: f64,
    pub weight_sum: f64,
    pub iterations: u64,
    pub running_means: Vec<Decimal>,
    pub match_accumulators: HashMap<u32, f64>,   // k matches → weighted sum
    pub per_individual: HashMap<String, f64>,    // id → weighted sum
}
```

`weighted_sum`, `weight_sum`, and the accumulator values use `f64` rather than `Decimal`. IS weights can exceed `Decimal`'s ~7.9×10²⁸ ceiling when a strong fixed bias is applied across many unknowns; `f64` avoids silent overflow. `running_means` stays as `Decimal` for the final precision-sensitive convergence check.

`running_mean() = weighted_sum / weight_sum` (importance-weighted mean).

### 9.2 Pedigree Probability Batch

`simulate_pedigree_probability_batch(pedigree, marker_set, params, bfs_layers, bias_value, batch_length, initial_sums, seed) -> BatchResult`

Estimates the average probability of the observed pedigree configuration (i.e., the likelihood of the known haplotypes given the pedigree structure).

```
Algorithm (per iteration):
1. Traverse unknowns in BFS order (root → leaves).
2. For each unknown U:
   a. Get parent's current haplotype.
   b. Compute biases for U using get_biases_for_individual().
   c. (new_hap, w, u) = mutate_haplotype(parent_hap, ...)
   d. Store new_hap in haplotypes map.
   e. w_total *= w       (product of biased transition probs)
   f. u_total *= u       (product of neutral transition probs)
3. Importance weight: iw = u_total / w_total
4. For each "non-simulated" edge (both endpoints have known haplotypes):
   p = get_edge_probability(parent, child, ...)
   If p == 0: ped_prob = 0 (this configuration is impossible)
5. ped_prob = Π(non-simulated edge probabilities)
6. Accumulate: (probability=ped_prob, importance_weight=iw)
```

**Importance weight derivation:**
```
The biased sampler draws haplotypes from q(h), not the neutral p(h).
The importance weight corrects: E_p[f] = E_q[f × p/q] = E_q[f × iw]
where iw = p(haplotype sequence) / q(haplotype sequence)
         = u_total / w_total
```

Batches run in parallel using Rayon's `par_iter`. Each parallel thread receives a distinct seed derived from `(trial_nr * 100 + model_index)`.

### 9.3 Matching Haplotypes Batch

`simulate_matching_haplotypes_batch(..., is_outside, ...) -> BatchResult`

The core computation: estimates the probability that the suspect haplotype originates from the pedigree.

```
Algorithm (per iteration):
1. Pick one "fixed" unknown U_f with probability proportional to picking_probabilities.
2. Assign suspect_haplotype to U_f in haplotypes map.
3. Traverse all OTHER unknowns in BFS order:
   a. Compute biases (treating U_f as Known).
   b. Mutate each unknown as before → accumulate w_total, u_total.
4. iw = u_total / w_total
5. Compute all edge probabilities → ped_prob = Π(non-simulated edges)
6. sim_prob_factor = Π(simulated edge probabilities)
   (only edges where the child was actually simulated in step 3)
7. fixed_prob = picking_probs[U_f] / Σ(picking_probs)
8. simulation_probability = fixed_prob × sim_prob_factor

For INSIDE calculation:
   non_excl_matches = count of non-excluded individuals whose haplotype
                      matches suspect_haplotype (including U_f)
   total_number     = total count including U_f (non-excluded)
   If non_excl_matches >= 1:
       cond = ped_prob / avg_pedigree_probability
       probability = cond / (simulation_probability × total_number)
       Accumulate: match_accumulators[non_excl_matches] += probability × iw
       For each non-excluded individual with matching haplotype:
           per_individual[id] += (probability / non_excl_matches) × iw
       weight_sum += iw

For OUTSIDE calculation:
   cond = ped_prob / avg_pedigree_probability
   probability = cond / simulation_probability
   Accumulate: weighted_sum += probability × iw; weight_sum += iw
```

**Mathematical foundation:**

The inside probability for exactly `k` matches is:
```
P(k matches | suspect) = E[ I(non_excl_matches = k) × ped_prob ]
                          / E[simulation_probability × k × avg_ped_prob]
```

The importance weighting (`iw = u/w`) ensures estimates are unbiased despite the biased sampler.

---

## 10. Three-Model Ensemble & Convergence

`simulation/convergence.rs` runs three independent simulation models with different bias levels and checks convergence of their grand mean.

### 10.1 Ensemble Trial Structure

```rust
pub struct EnsembleTrial {
    pub trial_nr: u32,
    pub model_results: [BatchResult; 3],
    pub converged: bool,
    pub grand_mean: Option<Decimal>,
}
```

`grand_mean() = mean(model_results[0].running_mean(), ..., model_results[2].running_mean())`

### 10.2 Convergence Check

```rust
fn check_convergence(results: &[BatchResult; 3], criterion: f64) -> bool
```

```
Algorithm:
1. grand_mean = mean(model[0].last_mean, model[1].last_mean, model[2].last_mean)
2. If grand_mean == 0: return false (never converge on zero)
3. For each model m:
   For each running_mean rm in model[m].running_means:
     If |rm - grand_mean| / grand_mean > criterion: return false
4. return true
```

All running means from all three models must stay within `criterion` fractionally of the grand mean. This is a strict criterion that requires all three models to have stabilised.

### 10.3 Ensemble Loop

```
Algorithm (run_ensemble_pedigree_probability / run_ensemble_matching_haplotypes):

criterion = params.convergence_criterion
tighten_count = 0

for trial_nr in 1.. (unbounded):
    If adaptive_bias: update model biases from schedule
    
    # Run 3 models in parallel (Rayon par_iter)
    for model in [0, 1, 2]:
        bias = per_model_bias[model]
        seed = trial_nr * 100 + model
        batch_result = simulate_*_batch(..., seed=seed)
    
    # Accumulate results cumulatively (carry across trials)
    for model in [0, 1, 2]:
        carry[model].weighted_sum += batch_result[model].weighted_sum
        carry[model].weight_sum   += batch_result[model].weight_sum
        carry[model].running_means.append(carry[model].running_mean())
    
    # Check for grand mean > 1 (pathological case)
    if grand_mean > 1 and tighten_count < 2:
        criterion /= 1.5
        tighten_count += 1
    
    # Check convergence
    if check_convergence(carry, criterion):
        if adaptive_bias: update schedule with current grand_means
        return EnsembleTrial { converged: true, grand_mean, ... }
    
    # Emit progress event (trial, model means, converged=false)
    if cancel_flag.load(): return Err(Cancelled)
```

Seeds are deterministic but differ across models and trials, ensuring the three parallel models explore different random paths.

---

## 11. Top-Level Simulation Flow

`simulation/mod.rs` orchestrates the entire pipeline.

```
run_simulation(pedigree, marker_set, params, progress_tx, cancel_flag):

Phase 0 — Classify individuals:
    If trace_mode:
        pedigree.add_trace(trace_haplotype)  → marks as Suspect
    Else:
        Find params.suspect in pedigree → mark as Suspect
    For each id in params.exclude: set individual.exclude = true

Phase 1 — Pedigree preparation:
    ext_ped = pedigree.clone()
    last_child_name = ext_ped.extend_pedigree()
    
    root_name     = pedigree.remove_irrelevant_individuals(inside=true, None)
    ext_root_name = ext_ped.remove_irrelevant_individuals(inside=false, Some(last_child_name))
    
    pedigree.reroot(root_name)
    ext_ped.reroot(ext_root_name)

Phase 2 — Picking probabilities:
    pedigree.calculate_picking_probabilities(trace_mode)
    ext_ped: only last_child_name gets probability = 1

Phase 3 — BFS layer setup:
    layers     = bfs_layers(pedigree, root_name)
    ext_layers = bfs_layers(ext_ped, ext_root_name)
    known_count = number of Known/Suspect individuals in pedigree

Step 1 — Average pedigree probability:
    If known_count <= 1:
        avg_ped_prob = Decimal::ONE  (trivially 1 for degenerate case)
    Else:
        trial = run_ensemble_pedigree_probability(pedigree, ...)
        avg_ped_prob = trial.grand_mean

Step 2 — Inside match probabilities (if not skip_inside):
    trial = run_ensemble_matching_haplotypes(pedigree, ..., is_outside=false)
    inside_probs = aggregate_match_counts(trial.model_results)
    per_individual = aggregate_per_individual(trial.model_results)

Step 3 — Outside match probability (if not skip_outside):
    ext_trial = run_ensemble_pedigree_probability(ext_ped, ...)
    ext_avg_ped_prob = ext_trial.grand_mean
    outside_trial = run_ensemble_matching_haplotypes(ext_ped, ..., is_outside=true)
    outside_prob = outside_trial.grand_mean

Return SimulationResult {
    inside_match_probabilities: { probabilities: {k→p}, average_pedigree_probability },
    outside_match_probability,
    per_individual_probabilities,
    trials, converged
}
```

---

## 12. File I/O

`matchy-io` handles all input/output with error type `IoError`.

### 12.1 Haplotype JSON Format

```json
{
  "IndividualName": {
    "DYS576": "15",
    "DYF387S1": "35;37.1"
  },
  "TRACE": { ... }
}
```

Allele string parsing:
```
"15"       → [(value=15, intermediate=None)]
"15.2"     → [(value=15, intermediate=Some(2))]
"15;16.3"  → [(value=15, intermediate=None), (value=16, intermediate=Some(3))]
```

If the key `"TRACE"` is present, the haplotype is stored as `pedigree.trace_haplotype` and returned separately. All matched individuals are promoted to `HaplotypeClass::Known`.

### 12.2 TGF (Trivial Graph Format)

```
<id> <name>
<id> <name>
...
#
<parent_id> <child_id>
<parent_id> <child_id>
...
```

Nodes before `#`, edges after. Simple and human-readable.

### 12.3 PLINK PED Format

Standard PLINK columns: `family_id individual_id paternal_id maternal_id sex phenotype`. Females (sex=2) are skipped. Only paternal relationships are imported.

### 12.4 Embedded Assets

Both `mutation_rates.csv` and `kits.json` are embedded into the binary at compile time using `include_str!`. This means the tool works standalone with no external data files.

**Kit format** (`kits.json`):
```json
{ "Y-37": ["DYS393", "DYS390", ...], "Y-67": [...] }
```

**Mutation rates** (`mutation_rates.csv`):
```
marker,mutation_rate
DYS393,0.00076
DYS390,0.00302
...
```

### 12.5 TOML / INI Config

The CLI accepts either TOML or legacy INI/CFG config files. The config is mapped directly to `SimulationParameters` via `From<Config>`.

---

## 13. Report Generation

`matchy-report/html.rs` renders HTML reports using embedded Minijinja templates.

Two report types:
- **Standard report** — for suspect mode, shows inside/outside probabilities, per-individual breakdown, convergence charts.
- **Trace report** — for trace mode, shows ranked per-individual probabilities, percentage relative to top match.

Key template data built in Rust:

```
haplotype_table_rows: per-marker rows across all individuals
    → highlight=true if any alleles differ (showing polymorphism)

marker_rows: sorted by name, showing rate and copy count

per_individual: sorted descending by probability

inside_k1_lr = (1 - f) / f  where f = inside_prob
    (likelihood ratio for at least one match)
```

The logo PNG is embedded as a base64 data URL (pure-Rust base64 encoder, no dependency):
```rust
fn base64_encode(data: &[u8]) -> String {
    // maps every 3 bytes to 4 base64 characters using the standard alphabet
}
```

---

## 14. CLI

`matchy-cli/main.rs` implements the command-line workflow:

```
1. Parse CLI args (config path, override flags).
2. load_config(config_path) → Config
3. Apply CLI flag overrides (skip_inside, skip_outside, trace_mode, adaptive_bias).
4. Resolve all file paths relative to config directory.
5. Load marker set: built-in kit name → load_kit() OR CSV file → read_marker_csv().
6. Load pedigree: .tgf → read_tgf(), .ped → read_ped().
7. Load haplotypes: read_haplotypes(pedigree, marker_set).
8. Config → SimulationParameters.
9. Spawn a thread to drain progress events (for capture into the HTML report charts).
10. run_simulation(pedigree, marker_set, params, Some(progress_tx), None).
11. print_summary(result) → stdout.
12. render_report() or render_trace_report() → HTML string (progress events passed for convergence charts).
13. Write {simulation_name}_report.html to results_path.
```

The CLI passes `Some(progress_tx)` for the progress channel so convergence chart data is captured for the HTML report. The cancel flag is `None` (no cancellation support in CLI mode).

---

## 15. Tauri Desktop App — Rust Backend

### 15.1 State

```rust
pub struct AppState {
    pub simulation: Mutex<Option<SimulationHandle>>,
}
pub struct SimulationHandle {
    pub cancel_flag: Arc<AtomicBool>,
}
```

`cancel_flag` is an `AtomicBool` shared between the handle (stored in `AppState`) and the simulation thread. The frontend calls `cancel_simulation` → Tauri → `handle.cancel()` → `flag.store(true)` → the simulation loop polls it and returns `Err(Cancelled)`.

### 15.2 Simulation Command

The Tauri `run_simulation` command runs the CPU-bound simulation on a blocking thread, emitting progress events via a Tauri event channel:

```
run_simulation command (async Tauri handler):
1. Cancel any existing simulation (via AppState).
2. Create new SimulationHandle, store in AppState, keep cancel_flag.
3. Parse inputs: TGF string → Pedigree, JSON string → haplotypes,
   kit name or CSV → MarkerSet.
4. Build SimulationParameters from request.
5. Spawn tokio::task::spawn_blocking for progress drainer:
       loop { if let Some(event) = rx.try_recv(): emit("simulation-progress", event) }
6. Spawn tokio::task::spawn_blocking for simulation:
       run_simulation(pedigree, marker_set, params, Some(tx), Some(cancel_flag))
7. Await simulation result.
8. On success: emit "simulation-complete", return SimulationResponse { success: true, result }.
9. On Cancelled: return { success: false, error: "cancelled" }.
```

### 15.3 Pedigree Commands

`build_extended_pedigree` runs the full prepare pipeline (extend, prune, reroot) and returns a `PedigreeData` JSON for the frontend to render as an SVG preview of the outside pedigree used in the simulation.

`validate_pedigree` checks both `is_dag` (no cycles) and `is_connected` (single connected component).

### 15.4 Report Command

`save_and_open_report` sanitises the filename to `[a-zA-Z0-9_-]+`, writes to the OS temp directory, then opens with the platform's default browser:
- Windows: `cmd /c start <path>`
- macOS: `open <path>`
- Linux: `xdg-open <path>`

### 15.5 Update Check Command

`check_for_updates(current_version)` queries `https://api.github.com/repos/genid/MatchY/releases/latest` via `reqwest` and returns:

```rust
pub struct UpdateInfo {
    pub latest_version: String,
    pub is_newer: bool,
    pub release_notes: String,
    pub download_url: String,
}
```

Semver comparison is a simple tuple `(major, minor, patch)`. Called silently on startup by the frontend; network errors are swallowed. If `is_newer` is true, the frontend shows an orange badge on the version number that opens a modal with release notes and a download link.

---

## 16. Frontend — React / TypeScript

### 16.1 Technology Stack

- **React 18** with TypeScript
- **React Router** — 5 top-level pages
- **Zustand** — global state store (no persistence; session files used instead)
- **@xyflow/react** (ReactFlow) + **dagre** — interactive DAG editor for pedigrees
- **Chart.js** + `chartjs-plugin-zoom` — live convergence charts
- **Tailwind CSS** — styling

### 16.2 Global State (`store/appStore.ts`)

```typescript
interface AppStore {
  pedigree: PedigreeData | null;
  pedigreeTgf: string;
  haplotypes: HaplotypesParseResult | null;
  haplotypesJson: string;
  selectedKitName: string | null;
  markers: MarkerInfo[];
  markerSetCsv: string | null;
  suspect: string | null;
  exclude: string[];
  simulation: {
    running: boolean;
    progress: ProgressEvent[];
    response: SimulationResponse | null;
    result: SimulationResult | null;
    error: string | null;
  };
}
```

### 16.3 `useSimulation` Hook

Manages the async lifecycle of a Tauri simulation invocation:

```typescript
async function startSimulation(params) {
  resetSimulation(); setSimulationRunning(true);
  const unlisten = await listen("simulation-progress", e => addProgressEvent(e.payload));
  try {
    const response = await invoke("run_simulation", { request });
    setSimulationResult(response);
  } catch (e) {
    if (!e.toString().includes("cancelled")) setSimulationError(e.toString());
  } finally {
    setSimulationRunning(false); unlisten();
  }
}
```

### 16.4 Pedigree Builder

The `PedigreeBuilder` page uses ReactFlow to provide an interactive DAG editor:

- **Nodes** are colored by `HaplotypeClass` (green=known, pink=suspect, blue=fixed, yellow=estimated, grey=unknown/excluded).
- **Connections** enforce single-parent: each node may have at most one parent.
- **Cycle detection** before accepting a connection:
  ```
  wouldCreateCycle(parentId, childId, relationships):
      BFS from parentId following parent-edges upward
      return childId is reachable (would create a cycle if edge is added)
  ```
- **Drag-to-empty-canvas**: dropping onto empty canvas opens a dialog to create and name a new node, automatically wiring the relationship.
- The raw TGF representation is kept in sync with each edit via `invoke("export_tgf")`.

### 16.5 Haplotype Editor

A spreadsheet-style table with:
- Rows = markers, Columns = individuals + optional TRACE column
- **Auto-save** (300ms debounce) on every cell change
- **Validation** per cell: must match `\d+(\.\d+)?(;\d+(\.\d+)?)*`; copy counts must be consistent across all individuals for each marker
- **Diversity highlighting**: markers where any allele values differ from another individual are highlighted amber
- **Keyboard navigation**: Tab/Shift+Tab, Enter/Shift+Enter for cell-to-cell movement
- **Per-individual CSV import**: auto-detects `,` vs `;` delimiter; maps `marker_name` / `alleles` column headers

### 16.6 Marker Sets Page

Provides two modes:
1. **Built-in kit** — select from embedded kits (RMplex, Yfiler plus, PowerPlex Y23, Combined)
2. **Custom set** — checkbox table of all known markers with per-marker rate and copy-count overrides

Auto-detects copy counts from loaded haplotypes:
```typescript
const detectedCopies = useMemo(() => {
  for each marker in haplotypes:
    count = max(semicolon-separated allele count across all individuals)
    if count > 1: detectedCopies[marker] = count
}, [haplotypes, allMarkers])
```

Displays overall mutation rate = `1 - Π(1 - mu_i)` for the selected set.

### 16.7 Convergence Charts

Each of the four simulation stages (pedigree probability, extended pedigree probability, inside match, outside match) gets a live-updating `ConvergenceChart` component:

- Filters `ProgressEvent[]` by stage
- Groups events by model index (0, 1, 2), showing three lines
- X-axis = batch number, Y-axis = running mean
- `animation: false` for smooth live updates
- Supports wheel/pinch zoom via `chartjs-plugin-zoom`
- `toBase64Image()` method (via `forwardRef`) used to capture chart images for the HTML report

### 16.8 SVG Pedigree Renderer (`utils/pedigreeSvg.ts`)

Renders the pedigree as a layered SVG for the HTML report:

```
Algorithm:
1. BFS from all roots to assign depth layers.
2. Within each layer, center nodes horizontally.
3. NODE_W=90, NODE_H=32, H_GAP=24, V_GAP=48, PADDING=20
4. Draw right-angle connector paths:
   M px,py L px,mid L cx,mid L cx,cy
   where mid = midpoint between parent and child layers
5. Color nodes by class:
   excluded → #e9ecf0 / dashed border
   known    → #c6f6d5 / green
   suspect  → #fed7d7 / pink
   unknown  → #e2e8f0 / dashed border
6. Encode as base64 data URL:
   "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(svg)))
```

### 16.9 Session Persistence (`utils/session.ts`)

Sessions bundle all pedigree, haplotype, marker set, simulation parameters, and results into a single `.matchy.json` file:

```typescript
interface SessionData {
  version: 1;
  pedigreeTgf: string;
  haplotypesJson: string;
  selectedKitName: string | null;
  markerSetCsv: string | null;
  suspect: string | null;
  exclude: string[];
  params: SessionParams;
  simulationResult: SimulationResult | null;
  simulationProgress: ProgressEvent[];
}
```

On load, all data is re-parsed through the Rust backend (TGF → Pedigree, JSON → haplotypes, etc.) to ensure consistency rather than deserialising potentially stale frontend state directly.
