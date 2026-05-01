![MatchY](python/assets/banner.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub Release](https://img.shields.io/github/v/release/genid/MatchY)](https://github.com/genid/MatchY/releases/latest)
[![GitHub issues](https://img.shields.io/github/issues/genid/MatchY)](https://github.com/genid/MatchY/issues)

MatchY is a forensic genetics tool for estimating match probabilities for Y-STR haplotypes. It uses Monte Carlo simulation with importance sampling, modelling mutations across pedigree structures to compute likelihood ratios for inside- and outside-pedigree hypotheses. Supports any number of markers, including multi-copy markers and intermediate alleles.

## Documentation

| Document | Description |
|----------|-------------|
| [GUI User Manual](docs/gui-manual.md) | Step-by-step guide for the desktop application |
| [CLI User Manual](docs/cli-manual.md) | Configuration reference, batch mode, and examples |
| [Technical Manual](docs/technical-manual.md) | Algorithm, statistical model, and Rust implementation |

## Download

Pre-built binaries are attached to every [GitHub Release](https://github.com/genid/MatchY/releases/latest).

### Desktop app

| Platform | File | Notes |
|----------|------|-------|
| Windows | `MatchY-{version}-x64-windows-setup.exe` | NSIS installer |
| Windows | `MatchY-{version}-x64-windows.msi` | MSI installer |
| Windows | `MatchY-{version}-x64-windows-portable.zip` | No install required |
| Linux | `MatchY-{version}-x86_64-linux.AppImage` | Make executable, then run |
| macOS ARM | `MatchY-{version}-aarch64-macos.dmg` | Unsigned — see note |

### CLI (zero dependencies)

| Platform | File |
|----------|------|
| Windows | `matchy-{version}-x86_64-windows.exe` |
| Linux | `matchy-{version}-x86_64-linux-musl` |
| macOS ARM | `matchy-{version}-aarch64-macos` |

> **macOS note:** Binaries are unsigned. Right-click → Open the first time to bypass Gatekeeper, or run `xattr -dr com.apple.quarantine <file>` in Terminal.

## CLI quick start

```bash
# Single analysis
matchy -c config.toml

# Batch: run all configs in a folder
matchy -c /path/to/cases/

# Trace mode (identify most likely donor)
matchy -c config.toml --trace-mode

matchy --help
```

See the [CLI manual](docs/cli-manual.md) for the full configuration reference.

## Build from source

### Desktop app (Windows)

```powershell
.\build.ps1
```

Requires: Rust stable, Node.js 20+. Artifacts appear in `C:\cargo-target\matchy\release\bundle\`.

### CLI

```bash
cargo build --release -p matchy-cli --manifest-path matchy/Cargo.toml
```

Static Linux binary (zero runtime dependencies):

```bash
rustup target add x86_64-unknown-linux-musl
cargo build --release -p matchy-cli --manifest-path matchy/Cargo.toml --target x86_64-unknown-linux-musl
```

## Python version (legacy)

A Python/Streamlit implementation is preserved in `python/` for reference. It is no longer the primary interface.

## License

MIT — see [LICENSE](LICENSE).

Copyright (c) 2026 Department of Pathology and Clinical Bioinformatics, Erasmus MC University Medical Center Rotterdam, The Netherlands; Institute of Medical Informatics and Statistics, Kiel University, University Hospital Schleswig-Holstein, Kiel, Germany; Chair of Epidemiology, Medical Biometry and Medical Informatics, Department of Medicine, Health and Medical University Erfurt, Erfurt, Germany.

## Citation

If you use MatchY in your research, please cite:

Zandstra D, Ralf A, Caliebe A, Nothnagel M, Krawczak M, Kayser M (2026). MatchY: A software implementation of pedigree-based calculation of Y-STR match probabilities. *Forensic Science International: Genetics*, 103518. https://doi.org/10.1016/j.fsigen.2026.103518
