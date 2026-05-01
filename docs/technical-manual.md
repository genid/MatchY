# MatchY — Technical Manual

## Overview

MatchY estimates Y-STR haplotype match probabilities by Monte Carlo simulation with importance sampling over a user-supplied pedigree. This document describes the statistical model, simulation algorithm, convergence criterion, and Rust implementation.

---

## Problem Statement

Given:
- A Y-chromosome pedigree (directed acyclic graph of individuals related by paternal descent)
- Known Y-STR haplotypes for a subset of individuals
- A Y-STR mutation model (per-marker rates, step-size distribution)

Compute:
1. **Pedigree probability** — P(all observed haplotypes | pedigree structure, mutation model)
2. **Inside-pedigree match probability** — P(at least one unknown non-excluded pedigree member matches the person of interest's haplotype)
3. **Outside-pedigree match probability** — P(a random individual just outside this pedigree matches the PoI haplotype)
4. **Likelihood ratios** derived from the above

---

## Mutation Model

### Single-copy markers

Each generation, an allele mutates with probability `μ` (the per-marker mutation rate). When a mutation occurs:

- With probability `1 − f_2`: the allele changes by ±1 repeat unit (one-step mutation).
- With probability `f_2`: the allele changes by ±2 repeat units (two-step mutation).

Direction (up or down) is equally likely for both step sizes.

The parameter `f_2` is the **two-step mutation fraction**, defaulting to 0.03.

### Multi-copy markers

For markers with `n` copies (e.g. DYF387S1 with 2 copies), each copy mutates independently with the same per-copy mutation rate. The observed haplotype is the multi-set of allele values across all copies.

### Matching

Two haplotypes match at a marker if their allele multi-sets are equal. A full haplotype match requires matching at every marker in the active set.

---

## Importance Sampling

Naive Monte Carlo estimates P(match) by simulating many random haplotypes for unknown individuals and counting how often the haplotype configuration is consistent with the observations. For rare events (low mutation rates, many generations), the naive estimator has very high variance — billions of samples would be needed.

MatchY uses **importance sampling** to reduce variance. Instead of sampling haplotypes uniformly, the simulation biases the sampling distribution toward configurations that are more likely to produce the observed outcome, then corrects for this bias using importance weights.

The bias parameter `b ∈ [0, 1)` controls the strength of the bias:
- `b = 0` → no bias (uniform sampling, equivalent to naive Monte Carlo)
- Higher `b` → stronger bias toward observed alleles, lower variance for rare events

The **adaptive bias** option (recommended) automatically tunes `b` during the simulation based on the effective sample size, finding an appropriate value without manual tuning.

### Importance weight

For each simulated haplotype configuration, the importance weight corrects for the sampling bias:

```
w = P(configuration | true mutation model) / P(configuration | biased model)
```

The estimator for probability `P` is:

```
P̂ = (1/N) Σ [ I(match) × w ]
```

where `I(match) = 1` if the simulated configuration produces a match and `0` otherwise.

---

## Three-Model Ensemble

To estimate convergence without a ground truth, MatchY runs three **independent** simulation models in parallel. Each model uses a different random seed and (slightly) different importance-sampling parameters.

The convergence criterion is met when the relative spread between the three models' current running means falls below the threshold `ε`:

```
spread = (max(P̂₀, P̂₁, P̂₂) − min(P̂₀, P̂₁, P̂₂)) / mean(P̂₀, P̂₁, P̂₂) ≤ ε
```

The default `ε = 0.02` (2%). The final reported probability is the mean of the three model estimates.

This approach is robust to local stochastic fluctuations: a single model might temporarily appear stable while still drifting, but three independent models drifting in the same direction simultaneously is rare.

---

## Pedigree Probability

The pedigree probability P(pedigree) is computed by summing over all possible haplotype assignments to unknown individuals, weighting each by its probability under the mutation model.

For a pedigree with root `r`, the probability factorises over edges:

```
P(all haplotypes) = P(root alleles) × ∏_{(parent,child) edges} P(child | parent)
```

where `P(child | parent)` is the probability of the child's haplotype given the parent's haplotype under the mutation model.

In importance sampling, the simulation draws haplotypes for unknown individuals from the biased distribution and accumulates weighted probability estimates.

---

## Inside-Pedigree Match Probability

The inside-pedigree match probability asks: given the PoI's haplotype, what is the probability that at least one other unknown, non-excluded individual in the pedigree has the same haplotype?

This is computed as:

```
P(inside match) = 1 − P(no other individual matches)
                = Σ_{k≥1} P(exactly k individuals match)
```

The simulation accumulates a count distribution `{k → P(k matches)}` and the reported probability is `Σ_{k≥1} P(k)`.

### Likelihood ratio (pedigree odds)

```
LR_pedigree = P(PoI is unique | observed) / P(at least one other matches | observed)
            = (1 − P_inside) / P_inside
```

This answers: how many times more probable is it that the PoI is the only matching individual versus at least one other pedigree member also matching?

---

## Outside-Pedigree Match Probability

The outside-pedigree match probability estimates P(random individual just outside this pedigree matches PoI), where "just outside" means: related to the pedigree but not typed or included.

The simulation extends the pedigree by adding one virtual individual as a child of each leaf node, then estimates the probability that this virtual individual matches the PoI haplotype. The result is averaged across all extended positions.

### Likelihood ratio (outside LR)

```
LR_outside = 1 / P_outside
```

This answers: how many times more probable is the evidence if the PoI is the trace donor than if a random outside individual is the donor?

---

## Trace Mode

In trace mode, there is no prior PoI — instead, the simulation identifies which pedigree member is the most likely donor of a crime-scene profile (the TRACE profile).

For each unknown individual `i`, the simulation estimates `P_i(TRACE match)` — the probability that individual `i` matches the TRACE haplotype. The results are reported as:
- Absolute match probabilities for each individual
- Relative ratios scaled to the most likely donor (100%)
- The outside-pedigree probability as a baseline (likelihood that the donor is not in the pedigree)

---

## Simulation Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| `batch_length` | N | 10,000 | Monte Carlo samples evaluated per batch. Each batch produces one data point on the convergence chart. |
| `convergence_criterion` | ε | 0.02 | Maximum relative spread between three model estimates. Simulation stops when this threshold is met. |
| `two_step_mutation_fraction` | f₂ | 0.03 | Fraction of mutations that are two-step (±2 repeat units). |
| `number_of_threads` | — | all CPUs | Rayon thread pool size. The three models run in parallel; within each model, per-batch parallelism uses the remaining threads. |
| `bias` | b | auto | Importance-sampling bias. Omit for automatic tuning. |
| `seed` | — | random | RNG seed for reproducibility. |

### Choosing convergence criterion

| ε | Precision | Use case |
|---|-----------|---------|
| 0.05 (5%) | Low | Exploratory, fast screening |
| 0.02 (2%) | Medium | Default, most casework |
| 0.01 (1%) | High | Publication, critical decisions |

Halving ε roughly quadruples the number of batches required.

### Choosing batch length

Larger batch lengths:
- Smoother convergence charts (less noise per data point)
- Fewer chart updates during a run
- No effect on final precision (precision is determined by ε)

The default of 10,000 is appropriate for most cases.

---

## Rust Implementation

### Crate structure

```
matchy/
├── crates/
│   ├── matchy-core/       # Core algorithm: simulation, mutation model, pedigree ops
│   ├── matchy-io/         # Config loading, pedigree/haplotype file formats, kits
│   ├── matchy-report/     # HTML report rendering (Minijinja templates)
│   ├── matchy-cli/        # CLI binary (clap argument parsing, batch runner)
│   └── matchy-app/        # Tauri desktop application + React/TypeScript frontend
└── Cargo.toml             # Workspace manifest (single version source)
```

### matchy-core

The core crate contains:

- **`simulation/mod.rs`** — Top-level `run_simulation()` entry point. Creates a Rayon thread pool sized to `number_of_threads`, then delegates to the convergence loop.
- **`simulation/convergence.rs`** — Convergence loop: runs batches until all three models meet the criterion. Emits `ProgressEvent` structs over an `mpsc` channel for live chart updates.
- **`simulation/importance.rs`** — Per-batch importance sampling. Splits `batch_length` samples across CPU threads (one RNG per thread to avoid per-iteration RNG initialisation overhead), accumulates weighted estimates.
- **`simulation/mutation.rs`** — Mutation model: generates child haplotypes from parent haplotypes under the single/two-step model.
- **`simulation/probability.rs`** — Pedigree probability accumulator.
- **`graph.rs`** — Petgraph-based DAG utilities: topological sort, BFS layers, ancestor/descendant queries, MRCA computation.
- **`models.rs`** — All public types: `SimulationParameters`, `SimulationResult`, `MatchProbabilities`, `StageStats`, `ProgressEvent`, etc.

### matchy-io

- **`config.rs`** — TOML and legacy INI config loader. Converts `Config` → `SimulationParameters`.
- **`tgf.rs`** / **`ped.rs`** — Pedigree file parsers.
- **`haplotypes_json.rs`** — JSON haplotype loader; assigns haplotypes to pedigree nodes and returns the optional TRACE profile.
- **`kits.rs`** — Loads embedded kit definitions from `assets/kits.json`.
- **`marker_csv.rs`** — Custom CSV marker set parser.
- **`mutation_rates.rs`** — Reads per-marker default mutation rates from `assets/mutation_rates.csv`.

### matchy-report

Uses **Minijinja** (a Rust Jinja2 port) to render HTML templates. Both `render_report()` and `render_trace_report()` accept a `SimulationResult` and optional JSON-serialised progress events, and return an HTML string. Templates are compiled into the binary via `include_str!`.

### matchy-cli

Thin binary that:
1. Parses CLI arguments with `clap`.
2. Detects single-file vs directory mode.
3. Calls `run_config()` for each config file.
4. Spawns a progress-collector thread to capture `ProgressEvent`s for report charts.
5. Writes the HTML report.

Uses **mimalloc** as the global allocator, which is critical for performance with the musl target — musl's default allocator serialises all allocations on a global mutex, causing severe slowdowns in parallel workloads.

### matchy-app

Tauri v2 desktop application. The Rust backend exposes Tauri commands (defined in `src/commands/`) that the React frontend calls via `invoke()`. Key commands:

| Command | Description |
|---------|-------------|
| `run_simulation` | Start a simulation; progress events streamed via Tauri events |
| `cancel_simulation` | Set the cancel flag |
| `generate_report` | Render and return HTML |
| `save_and_open_report` | Write HTML to disk and open in browser |
| `parse_tgf` / `parse_ped` | Parse uploaded pedigree files |
| `export_tgf` | Serialize current pedigree to TGF |
| `validate_pedigree` | Check DAG, connectivity |
| `build_extended_pedigree` | Extend pedigree for outside calculation |
| `parse_haplotypes_json` / `export_haplotypes_json` | Haplotype serialisation |
| `list_kits` / `load_kit` / `load_custom_csv` | Marker set management |

### Parallelism

All parallelism uses **Rayon**. The thread pool is created explicitly with `ThreadPoolBuilder` so its size is controlled by the user parameter, not Rayon's global default.

Within each batch:
```
n_chunks = rayon::current_num_threads()
chunk_size = batch_length / n_chunks
```
Each chunk runs sequentially with a single RNG (seeded once per chunk). This avoids the cost of initialising a ChaCha20 RNG state for every iteration — an initialisation that was the root cause of a 10–16× performance regression when per-iteration RNGs were used.

The three convergence models themselves also run in parallel, so the effective thread utilisation is:
```
3 models × (batch_length / chunk_size) chunks each
```

### Reproducibility

When a `seed` is specified, the simulation is fully deterministic: the same seed with the same parameters on the same binary produces bit-identical results. Without a seed, a fresh random seed is generated at startup.

---

## Performance Notes

| Factor | Impact | Recommendation |
|--------|--------|---------------|
| Thread count | Linear speedup up to ~16 threads for typical pedigrees | Use all available cores |
| Batch length | No effect on wall time per batch; larger = less chart overhead | 10,000 is optimal for most cases |
| Convergence criterion | Quadratic: halving ε ≈ 4× more batches | Use 0.02 for casework |
| Allocator (Linux/musl) | musl's default allocator serialises parallel allocations | Use provided musl binary (mimalloc) |
| Pedigree size | Larger pedigrees = more unknowns to simulate = slower per batch | No workaround; unavoidable |
| Marker count | Linear: more markers = more per-marker work | Use targeted panels for speed |

---

## Allele Format Reference

| Format | Example | Meaning |
|--------|---------|---------|
| Integer | `14` | Standard STR repeat count |
| Decimal | `9.3` | Intermediate allele (e.g. 9 repeats + partial) |
| Multiple values | `["35", "37"]` | Multi-copy marker (two distinct copies) |

Intermediate alleles are treated as distinct allele values — a mutation from 9.3 goes to 8.3 or 10.3 (one-step), not to 9 or 10.

---

## Reported Statistics

### StageStats (per simulation stage)

| Field | Description |
|-------|-------------|
| `iterations_per_model` | Total MC samples used by each of the three models |
| `model_probabilities` | Final probability estimate from each model |
| `total_iterations()` | Sum across all models |
| `formatted_runtime()` | Wall-clock time for this stage (s / m s / h m / d h) |

### SimulationResult

| Field | Description |
|-------|-------------|
| `inside_match_probabilities` | `MatchProbabilities` struct with count distribution and average pedigree probability |
| `outside_match_probability` | Scalar `Decimal` |
| `per_individual_probabilities` | Map of individual ID → match probability (trace mode only) |
| `trials` | Number of batches until convergence |
| `converged` | Whether the convergence criterion was met |
| `total_runtime_secs` | Total wall-clock time in seconds |
