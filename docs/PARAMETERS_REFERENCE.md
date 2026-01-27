# MatchY Parameters Reference Guide

Complete reference of all hard-coded constants, default parameter values, flags, and configuration defaults.

---

## 1. SIMULATION PARAMETERS

### Core Simulation Settings

| Parameter | Default Value | Location | Description | Effect of Increasing | Effect of Decreasing |
|-----------|--------------|----------|-------------|---------------------|---------------------|
| **two_step_mutation_fraction** | `0.03` | config.ini:14, Home.py:149-158 | Probability modifier for 2-step mutations relative to 1-step | Higher values = more 2-step mutations allowed, more permissive matching | Lower values = fewer 2-step mutations, stricter matching |
| **batch_length** | `10000` | config.ini:15, Home.py:160-166 | Number of iterations per trial for model convergence | More iterations = better convergence, longer runtime | Fewer iterations = faster but less reliable results |
| **convergence_criterion** | `0.02` | config.ini:16, Home.py:168-175 | Maximum relative difference (2%) between model means for validity | Higher threshold = models converge easier, less strict | Lower threshold = models must agree more closely, stricter convergence |
| **bias** | `None` | config.ini:17, SimulationParameters:955 | Importance sampling bias value (None = auto-calculate) | Higher bias = stronger direction preference in mutations | Lower bias = weaker preference, closer to uniform sampling |
| **number_of_threads** | `1` (functions), `min(4, CPU_count)` (GUI) | Home.py:237-239 | Number of parallel threads for simulation | More threads = faster execution (up to CPU limit) | Fewer threads = slower but less resource intensive |

**Valid Ranges:**
- `two_step_mutation_fraction`: [0.0, 1.0]
- `batch_length`: [0, ∞) (practically should be ≥ 1000)
- `convergence_criterion`: [0.0, ∞) (typical: 0.001-0.1)
- `bias`: [0.0, 0.5] or None
- `number_of_threads`: [1, CPU_count]

---

## 2. ADAPTIVE BIAS SYSTEM

### Adaptive Bias Schedule

| Trial/Rank | Bias Value | Location | Description | Effect of Changing |
|------------|-----------|----------|-------------|-------------------|
| **Trial 1 (all models)** | `0.10` | simulation.py:243 | Initial bias for all 3 models in first trial | Higher = stronger initial preference, lower = more exploration |
| **Best Model** | `0.05` | simulation.py:258 | Bias for best performing model | Lower = less guidance (more exploration), higher = more guidance |
| **Middle Model** | `0.15` | simulation.py:258 | Bias for middle performing model | Affects balance between exploration/exploitation |
| **Worst Model** | `0.25` | simulation.py:258 | Bias for worst performing model | Higher = strong compensation for poor performance |

**How Adaptive Bias Works:**
- Ranks 3 models by performance each trial
- Best model gets lowest bias (explores less, exploits current good strategy)
- Worst model gets highest bias (explores more to find better strategy)
- Only active when `adaptive_bias=True` flag is set

---

## 3. DEFAULT BIAS CALCULATION (When bias=None)

| Component | Formula | Location | Description |
|-----------|---------|----------|-------------|
| **Base Formula** | `min(max(0.1, (0.8 / (1 + distance_to_mrca))), 0.4)` | models.py:913 | Auto-calculated bias based on pedigree distance |
| **Minimum Bias** | `0.1` | models.py:913 | Lower bound (10%) |
| **Maximum Bias** | `0.4` | models.py:913 | Upper bound (40%) |
| **Distance Factor** | `0.8 / (1 + distance)` | models.py:913 | Inverse relationship with MRCA distance |
| **Disable Bias** | `bias_value <= 0` | models.py:914-915 | Returns empty bias list (no biasing) |

**Effect of Distance to MRCA:**
- Distance 0: bias = 0.4 (clamped to max)
- Distance 1: bias = 0.4 (0.8/2 = 0.4)
- Distance 3: bias = 0.2 (0.8/4 = 0.2)
- Distance 7: bias = 0.1 (0.8/8 = 0.1, clamped to min)
- Distance 20+: bias = 0.1 (clamped to min)

---

## 4. BATCH PROCESSING PARAMETERS

| Parameter | Value | Location | Description | Effect of Changing |
|-----------|-------|----------|-------------|-------------------|
| **Min Batch Size** | `256` | simulation.py:236 | Minimum iterations per batch | Larger = fewer context switches, more memory | Smaller = more frequent updates, less memory |
| **Max Batch Size** | `512` | simulation.py:236 | Maximum iterations per batch | Larger = less overhead, longer between updates | Smaller = more frequent progress updates |
| **Batch Calculation** | `window // (threads * 4)` | simulation.py:236 | Formula for dynamic batch size | Automatically scales with window size and threads |
| **Write Buffer** | `100` | simulation.py:313-315 | Iterations between disk writes | Larger = less I/O overhead, larger memory buffer | Smaller = more frequent saves, less data loss risk |

**Batch Size Formula:**
```python
batch_size = min(max(256, window // max(1, threads * 4)), 512)
```
- With 10,000 window & 4 threads: batch_size = 512 (capped at max)
- With 2,000 window & 4 threads: batch_size = 312
- With 500 window & 4 threads: batch_size = 256 (capped at min)

---

## 5. MODEL CONVERGENCE & RETRY LOGIC

| Parameter | Value | Location | Description | Effect of Changing |
|-----------|-------|----------|-------------|-------------------|
| **Number of Models** | `3` | simulation.py:227-232 | Ensemble model count (HARDCODED) | Cannot change without major code refactoring |
| **Max Tightening Count** | `2` | simulation.py:333 | Maximum retry attempts with stricter thresholds | More retries = more chances but longer runtime |
| **First Retry Divisor** | `1.5` | simulation.py:337 | Threshold reduction on first retry | Larger = gentler relaxation, smaller = more aggressive |
| **Second Retry Divisor** | `3.375` | simulation.py:926 | Threshold reduction on second retry | Larger = gentler relaxation, smaller = more aggressive |

**Retry Behavior:**
- Attempt 1: Use `convergence_criterion` as configured
- Attempt 2: Use `convergence_criterion / 1.5`
- Attempt 3: Use `convergence_criterion / 1.5 / 1.5` = threshold / 2.25
- Outside calc retry: Use `convergence_criterion / 3.375`

---

## 6. BOOLEAN FLAGS & MODES

| Flag | Default | Location | Description | True Effect | False Effect |
|------|---------|----------|-------------|------------|--------------|
| **trace_mode** | `False` | main.py:171, Home.py:70 | Enable trace donor identification mode | Analyzes trace vs all pedigree members | Analyzes specific suspect |
| **skip_inside** | `False` | main.py:154, Home.py:177 | Skip inside pedigree probability | Skips P(match\|inside pedigree) calculation | Calculates inside probability (normal) |
| **skip_outside** | `False` | main.py:162, Home.py:183 | Skip outside pedigree probability | Skips P(match\|outside pedigree) calculation | Calculates outside probability (normal) |
| **adaptive_bias** | `False` | main.py:180, Home.py:189 | Enable adaptive bias adjustment | Uses dynamic per-model bias (0.05/0.15/0.25) | Uses fixed bias from config |
| **exclude** | `False` | models.py:173 | Mark individual as excluded from analysis | Individual skipped in calculations | Individual included normally |
| **is_outside** | `False` | simulation.py:220 | Flag for outside probability calculation | Uses extended pedigree | Uses original pedigree |

---

## 7. GUI INPUT CONSTRAINTS

### Number Input Ranges

| Input Field | Min | Max | Step | Default | Location |
|-------------|-----|-----|------|---------|----------|
| **Two-Step Mutation Fraction** | 0.0 | 1.0 | 0.01 | 0.03 | Home.py:149-158 |
| **Batch Length** | 0 | ∞ | 1 | 10000 | Home.py:160-166 |
| **Convergence Criterion** | 0.0 | ∞ | 0.0001 | 0.02 | Home.py:168-175 |
| **Bias Value** | 0.0 | 1.0 | 0.01 | None | Home.py:206-215 |
| **Number of Threads** | 1 | CPU_count | 1 | min(4, CPU) | Home.py:237-239 |
| **Marker Mutation Rate** | 0.0000 | 1.0000 | 0.0001 | 0.0001 | Marker_sets.py:57-60 |

---

## 8. PEDIGREE VISUALIZATION CONSTANTS

| Constant | Value | Location | Description | Usage |
|----------|-------|----------|-------------|-------|
| **Window Width** | `1000` | config.ini:2 | Pedigree diagram width (pixels) | GUI display size |
| **Window Height** | `700` | config.ini:3 | Pedigree diagram height (pixels) | GUI display size |
| **Known Node Color** | `#b2d3c2` | config.ini:5 | Light green for known haplotypes | Visual distinction |
| **Unknown Node Color** | `#eeeeee` | config.ini:6 | Light gray for unknown haplotypes | Visual distinction |
| **Suspect Node Color** | `#ff0000` | config.ini:7 | Red for suspect individual | Highlight target |
| **Excluded Unknown Color** | `#888888` | config.ini:8 | Dark gray for excluded unknowns | Visual distinction |
| **Excluded Known Color** | `#5a7d6c` | config.ini:9 | Dark green for excluded knowns | Visual distinction |
| **Edge Color** | `#aaaaaa` | config.ini:10 | Gray for relationship lines | Neutral connections |

---

## 9. MUTATION PROBABILITY DISTRIBUTION

| Mutation Step | Probability | Location | Description |
|---------------|------------|----------|-------------|
| **No Mutation (0)** | `1.0 - mu` | simulation.py:48 | Probability of no change |
| **1-Step Up (+1)** | `(mu - two_step) / 2.0` | simulation.py:49 | Half of remaining probability after 2-step |
| **1-Step Down (-1)** | `(mu - two_step) / 2.0` | simulation.py:49 | Symmetric with +1 |
| **2-Step Up (+2)** | `(mu * two_step_factor) / 2.0` | simulation.py:49 | Half of 2-step probability |
| **2-Step Down (-2)** | `(mu * two_step_factor) / 2.0` | simulation.py:49 | Symmetric with +2 |

**Symmetric Distribution:** Division by `2.0` ensures equal probability for up/down mutations.

**Example with mu=0.004, two_step_factor=0.03:**
- 2-step total: 0.004 × 0.03 = 0.00012
- 1-step total: 0.004 - 0.00012 = 0.00388
- P(no mutation) = 0.996
- P(+1 or -1) = 0.00194 each
- P(+2 or -2) = 0.00006 each

---

## 10. EXTENDED PEDIGREE PARAMETERS

| Parameter | Value | Location | Description | Purpose |
|-----------|-------|----------|-------------|---------|
| **Last Child Picking Prob** | `Decimal(1)` | simulation.py:764 | 100% probability for last child in extended pedigree | Forces outside match calculation through specific path |
| **Other Individuals Prob** | `Decimal(0)` | simulation.py:769 | 0% probability for other extended individuals | Ensures calculations focus on target path |

---

## 11. DATA CLASS DEFAULTS

| Class | Field | Default | Location | Description |
|-------|-------|---------|----------|-------------|
| **Individual** | `haplotype` | `Haplotype()` | models.py:171 | Empty haplotype object |
| **Individual** | `exclude` | `False` | models.py:173 | Not excluded by default |
| **Individual** | `picking_probability` | `None` | models.py:175 | Auto-calculated if needed |
| **Individual** | `closest_known_individuals` | `[]` | models.py:176 | Empty list |
| **Individual** | `closest_known_distance` | `None` | models.py:177 | Not calculated initially |
| **Pedigree** | `individuals` | `[]` | models.py:255 | Empty list |
| **Pedigree** | `relationships` | `[]` | models.py:256 | Empty list |
| **Pedigree** | `picking_probabilities` | `{}` | models.py:257 | Empty dict |
| **SimulationParameters** | `bias` | `None` | models.py:955 | Auto-calculate |
| **SimulationParameters** | `user_name` | `None` | models.py:956 | Optional |
| **Marker** | `number_of_copies` | `None` | models.py:50 | Optional |
| **Allele** | `intermediate_value` | `None` | models.py:82 | Optional (for .1 notation) |
| **Allele** | `mutation_value` | `None` | models.py:83 | Optional |
| **Allele** | `mutation_probability` | `None` | models.py:84 | Optional |

---

## 12. CLI & FILE DEFAULTS

| Parameter | Default Value | Location | Description |
|-----------|--------------|----------|-------------|
| **Config File Path** | `"config.ini"` | main.py:145 | Default config file name |
| **Simulation Name** | `""` (empty) | Home.py:248 | No default name |
| **User Name** | `""` (empty) | Home.py:249 | No default user |
| **Suspect Index** | `0` | Home.py:91 | First individual in list |
| **Excluded List** | `[]` | Home.py:109 | No exclusions |
| **Haplotype Class** | `"unknown"` | models.py:172 | Default individual type |
| **File Uploader Key** | `0` | Haplotypes_editor.py:112 | Initial counter |

---

## 13. UI DISPLAY CONSTANTS

| Constant | Value | Location | Description |
|----------|-------|----------|-------------|
| **CLI Output Height** | `400` px | Home.py:410 | Text area for CLI output |
| **Log Container Height** | `600` px | Home.py:511 | Progress/log display |
| **Container Border** | `True` | Home.py:511 | Enable border on containers |

---

## 14. STRING CONSTANTS & ENUMS

### Haplotype Class Values
- `"unknown"` - Individual with unknown haplotype (default)
- `"known"` - Individual with known haplotype
- `"suspect"` - Individual marked as suspect/trace

### Bias Direction Values
- `"up"` - Mutation increases allele value
- `"down"` - Mutation decreases allele value
- `"none"` - No mutation occurred

### Execution Modes (Home.py:225-227)
- `"Run in Browser (Single-threaded)"` - Default Streamlit execution
- `"Submit to CLI (Multi-threaded)"` - Subprocess execution with multiprocessing

---

## 15. VALIDATION RULES

| Validation | Rule | Location | Error Behavior |
|------------|------|----------|----------------|
| **Bias Range** | `0.0 <= bias <= 0.5` | config.py:45-46 | Raises ValueError |
| **Thread Count** | `1 <= threads <= CPU_count` | Home.py:237-238 | Clamped to range |
| **Two-Step Fraction** | `0.0 <= factor <= 1.0` | Home.py:151-152 | Clamped to range |
| **Convergence Criterion** | `threshold >= 0.0` | Home.py:170 | Clamped to min |
| **Batch Length** | `window >= 0` | Home.py:162 | Clamped to min |

---

## 16. CRITICAL HARDCODED VALUES (Cannot Change Without Code Refactoring)

| Constant | Value | Impact | Location |
|----------|-------|--------|----------|
| **Number of Models** | `3` | Ensemble averaging uses exactly 3 models | simulation.py:227-232 |
| **Probability Distribution Divisor** | `2.0` | Ensures symmetric up/down mutations | simulation.py:49 |
| **Decimal Precision** | Python `Decimal` class | High-precision probability calculations | Throughout |
| **Model Ranking Count** | `3` positions | Best/Middle/Worst ranking | simulation.py:256-258 |

---

## QUICK REFERENCE: Most Important Parameters

| Parameter | Critical For | Typical Range | Recommendation |
|-----------|--------------|---------------|----------------|
| **stability_window** (Batch length) | Result reliability | 5,000-20,000 | Use 10,000 for standard cases |
| **model_validity_threshold** (Convergence criterion) | Convergence strictness | 0.01-0.05 | Use 0.02 for balance |
| **two_step_mutation_factor** (Two-step mutation fraction) | Mutation model | 0.01-0.05 | Use 0.03 (literature standard) |
| **adaptive_bias** | Performance | True/False | Use True for complex pedigrees |
| **bias** | Sampling efficiency | None or 0.1-0.3 | Use None for auto-calculation |
| **number_of_threads** | Runtime | 1-CPU_count | Use min(4, CPU) for GUI |

---

## PERFORMANCE TUNING GUIDE

### For Faster Results (Less Accurate)
- Decrease `stability_window` (Batch length) to 5000
- Increase `model_validity_threshold` (Convergence criterion) to 0.05
- Disable `adaptive_bias`
- Use fewer `number_of_threads` (reduces overhead for small jobs)

### For More Accurate Results (Slower)
- Increase `stability_window` (Batch length) to 20000+
- Decrease `model_validity_threshold` (Convergence criterion) to 0.01
- Enable `adaptive_bias`
- Increase `number_of_threads` to max CPU count

### For Complex Pedigrees
- Enable `adaptive_bias` = True
- Increase `stability_window` (Batch length) to 15000+
- Use lower `model_validity_threshold` (Convergence criterion) (0.015)
- Consider custom `bias` values if auto-calculation fails

---

*Last Updated: 2025-12-19*
*MatchY Version: trace_profiles branch*
