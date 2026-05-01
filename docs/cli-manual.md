# MatchY — CLI User Manual

## Overview

The MatchY CLI (`matchy`) is a self-contained command-line binary for running Y-STR match probability simulations. It produces the same HTML reports as the desktop application and is suitable for batch processing, server environments, and automated pipelines.

The Linux binary is statically linked (musl) and runs on any x86-64 Linux system without installing dependencies.

---

## Installation

Download the appropriate binary from the [GitHub Releases](https://github.com/genid/MatchY/releases/latest) page:

| Platform | File |
|----------|------|
| Windows | `matchy-{version}-x86_64-windows.exe` |
| Linux (any distro) | `matchy-{version}-x86_64-linux-musl` |
| macOS Intel | `matchy-{version}-x86_64-macos` |
| macOS ARM | `matchy-{version}-aarch64-macos` |

Make the binary executable on Linux/macOS:

```bash
chmod +x matchy-{version}-x86_64-linux-musl
mv matchy-{version}-x86_64-linux-musl /usr/local/bin/matchy
```

> **macOS:** Binaries are unsigned. Run `xattr -dr com.apple.quarantine matchy` to bypass Gatekeeper.

---

## Quick Start

```bash
# Run a single simulation
matchy -c config.toml

# Run all configs in a folder (batch mode)
matchy -c /path/to/cases/

# Trace mode (identify most likely donor)
matchy -c config.toml --trace-mode

# Skip outside-pedigree stage
matchy -c config.toml --skip-outside

# Show help
matchy --help
```

---

## Command-Line Arguments

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--config` | `-c` | `config.toml` | Path to a config file (`.toml` or `.ini`) **or** a folder of config files. |
| `--skip-inside` | `-i` | off | Skip the within-pedigree match probability stage. |
| `--skip-outside` | `-o` | off | Skip the outside-pedigree random match probability stage. |
| `--trace-mode` | `-t` | off | Identify the most likely trace donor from pedigree members. |
| `--adaptive-bias` | `-a` | off | Enable adaptive importance-sampling bias tuning. |
| `--help` | `-h` | — | Show help message. |
| `--version` | `-V` | — | Show version. |

CLI flags **override** the corresponding settings in the config file. This lets you use a shared base config and apply per-run overrides without editing the file.

---

## Configuration File

MatchY uses TOML as its primary config format. Legacy `.ini` files from the Python version are also supported and detected automatically.

### Full TOML reference

```toml
[simulation]
name = "my-analysis"          # Simulation name (used in report filename)
user = "Jane Smith"           # Analyst name (shown in report header)
results_path = "./results"    # Output directory for the HTML report

[files]
pedigree = "pedigree.tgf"              # Path to pedigree file (.tgf or .ped)
known_haplotypes = "haplotypes.json"   # Path to haplotypes JSON

# Use ONE of the following:
marker_set = "RMplex"                  # Built-in kit name (see list below)
# marker_set_file = "markers.csv"      # Custom CSV marker set

[parameters]
two_step_mutation_fraction = 0.03   # Fraction of two-step mutations (default: 0.03)
batch_length = 10000                # MC samples per batch (default: 10000)
convergence_criterion = 0.02        # Max inter-model relative spread (default: 0.02)
number_of_threads = 8               # CPU threads (default: all available)
bias = 0.5                          # Fixed IS bias (omit for automatic)
seed = 42                           # Random seed (omit to use default seed 0)

[mode]
suspect = "John"          # Person of interest ID (must match pedigree node name)
exclude = ["Uncle", "B"]  # Individuals excluded from match calculation
skip_inside = false       # Skip within-pedigree stage
skip_outside = false      # Skip outside-pedigree stage
trace_mode = false        # Identify most likely donor (requires TRACE in JSON)
adaptive_bias = false     # Differentiate IS bias across ensemble models (default: off)
```

### Minimal TOML (required fields only)

```toml
[simulation]
name = "my-analysis"

[files]
pedigree = "pedigree.tgf"
known_haplotypes = "haplotypes.json"
marker_set = "RMplex"

[mode]
suspect = "John"
```

### Embedded kits

| Kit name | Description |
|----------|-------------|
| `RMplex` | 30-marker Rapidly Mutating (RM) Y-STR panel |
| `Yfiler plus` | Promega Yfiler Plus — 27 Y-STR markers |
| `PowerPlex Y23` | Promega PowerPlex Y23 — 23 Y-STR markers |
| `Combined` | Union of RMplex, Yfiler Plus, and PowerPlex Y23 |

### Legacy INI format

The original Python `.ini` format is still supported:

```ini
[pedigree]
simulation_name = my-analysis
path = pedigree.tgf
known_haplotypes = haplotypes.json
marker_set = RMplex
suspect = John
two_step_mutation_fraction = 0.03
batch_length = 10000
convergence_criterion = 0.02
```

---

## Input File Formats

### Pedigree (`.tgf`)

Tab-separated graph format. Nodes are listed before the `#` separator, edges after.

Node format: `ID<tab>Name` — integer IDs are conventional and required by some external TGF tools. Edges reference nodes by ID.

```
1	Father
2	Son1
3	Suspect
#
1	2
1	3
```

The parser also accepts names directly as IDs (as exported by the GUI), but integer IDs are recommended for hand-authored files.

### Pedigree (`.ped`)

Standard 6-column PED format: `FamilyID IndividualID FatherID MotherID Sex Phenotype`. The Y-chromosome lineage is inferred from paternal IDs.

### Haplotypes (`.json`)

```json
{
  "Father": {
    "DYS19": ["14"],
    "DYS389I": ["13"],
    "DYS389II": ["29"],
    "DYF387S1": ["35", "37"]
  },
  "Suspect": {
    "DYS19": ["14"],
    "DYS389I": ["13"]
  },
  "TRACE": {
    "DYS19": ["14"],
    "DYS389I": ["13"]
  }
}
```

- Individual IDs must match the pedigree node names.
- Marker names must match those in the active marker set.
- Multi-copy markers have multiple values in the array.
- The special key `"TRACE"` is used as the crime-scene profile in trace mode.
- Individuals with no entry are treated as unknown (haplotype to be simulated).

### Custom marker set (`.csv`)

```csv
marker_name,mutation_rate,copies
DYS19,0.002,1
DYS389I,0.003,1
DYF387S1,0.004,2
```

- `copies` is optional; if omitted it is inferred from the haplotype data.
- Mutation rates are per-generation, single-step.

---

## Output

The CLI writes an HTML report to `results_path/{simulation_name}_report.html`. The report is identical to those produced by the desktop application and includes:

- Match probability summary table
- Narrative interpretation
- Convergence charts
- Performance metrics
- Simulation parameters
- Pedigree diagram
- Known haplotype tables

A summary is also printed to stdout:

```
=== MatchY Simulation Results ===
Simulation : my-analysis
Converged  : true (trials: 47)

Inside-pedigree match probabilities:
  Average pedigree probability: 1.23E-5
  P(at least 1 match) = 0.023456

Outside-pedigree match probability:
  P(match outside) = 0.00142

Per-individual match probabilities:
  Son1     : 0.01823  (LR = 54.8)
  Uncle    : 0.00412  (LR = 243)
```

**Per-individual match probability** is P(individual has the same haplotype as the PoI), estimated by Monte Carlo simulation. The **LR** for each individual is `1 / P(match)` and is the primary forensic metric: it quantifies how many times more likely an observed match is under the hypothesis that this individual is the true donor, compared to coincidental haplotype sharing.

The **pedigree odds** shown in the report is an odds ratio — P(PoI is the only match) / P(at least one other also matches) — and should not be confused with a likelihood ratio.

---

## Batch Mode

Pass a directory to `-c` to run all `.toml` and `.ini` files it contains (`.cfg` is also accepted as an alias for `.ini`), in alphabetical order:

```bash
matchy -c /path/to/cases/
```

Output:

```
INFO  Batch mode: 5 config(s) in "/path/to/cases/"
INFO  [1/5] Running "case_01.toml"
...
INFO  [5/5] Running "case_05.toml"
INFO  Batch complete: all 5 succeeded.
```

If one config fails (missing file, bad format, etc.), the error is logged and the remaining configs continue:

```
ERROR Failed: Cannot open pedigree "missing.tgf"
WARN  Batch complete: 1/5 failed.
```

CLI flags apply to all configs in the batch:

```bash
# Run all cases in trace mode, skipping outside stage
matchy -c /path/to/cases/ --trace-mode --skip-outside
```

---

## Log Output

The CLI uses structured logging via `tracing`. The default level is `INFO`. To change:

```bash
RUST_LOG=debug matchy -c config.toml    # verbose
RUST_LOG=warn  matchy -c config.toml    # warnings and errors only
RUST_LOG=error matchy -c config.toml    # errors only
```

Key log messages:

| Message | Meaning |
|---------|---------|
| `Simulation 'X' — N thread(s), batch_length=Y, convergence_criterion=Z` | Startup summary |
| `Pedigree probability converged after N trials: X` | Pedigree stage done |
| `Match probability converged after N trials: X` | Inside/outside stage done |
| `Report written to "..."` | Output file location |

---

## Examples

The `examples/` directory contains ready-to-run configurations:

```bash
# Standard suspect analysis
matchy -c examples/mockfo.toml

# Trace mode (identify most likely donor)
matchy -c examples/mockfo_trace.toml

# Custom marker set
matchy -c examples/simple_suspect.toml
```

---

## Building from Source

```bash
# Standard build
cargo build --release -p matchy-cli --manifest-path matchy/Cargo.toml

# Static Linux binary (no runtime dependencies)
rustup target add x86_64-unknown-linux-musl
cargo build --release -p matchy-cli --manifest-path matchy/Cargo.toml \
    --target x86_64-unknown-linux-musl
```
