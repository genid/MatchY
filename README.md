# MatchY

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-available-blue.svg)](https://hub.docker.com/r/dionzand/matchy)
[![Docker CLI](https://img.shields.io/badge/docker%20cli-available-blue.svg)](https://hub.docker.com/r/dionzand/match-cli)
[![GitHub issues](https://img.shields.io/github/issues/genid/MatchY)](https://github.com/genid/MatchY/issues)
[![GitHub stars](https://img.shields.io/github/stars/genid/MatchY)](https://github.com/genid/MatchY/stargazers)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](CONTRIBUTING.md)

MatchY is a powerful pedigree-based tool designed to estimate match probabilities for Y-STR haplotypes. Its mathematical framework leverages marker mutation rates, pedigree structure, and the known haplotypes of individuals within the family tree. By combining this data, the tool accurately estimates match probabilities with a person of interest using a Monte Carlo simulation with importance sampling to model mutations.

MatchY supports any number of Y-STR markers, including multi-copy markers and intermediate alleles.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
  - [Docker images (preferred)](#docker-images-preferred)
  - [Local installation (Linux only)](#local-installation-linux-only)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)
- [Support](#support)

## Features

- **Monte Carlo Simulation with Importance Sampling**: Accurately models Y-STR mutations across generations
- **Flexible Marker Support**: Works with any number of Y-STR markers, including multi-copy markers and intermediate alleles
- **Pedigree-Based Analysis**: Leverages family tree structure for precise probability estimates
- **Multiple Interfaces**:
  - Interactive web-based GUI (Streamlit)
  - Command-line interface for batch processing
- **Comprehensive Reporting**: Generate detailed HTML/PDF reports with visualizations
- **Trace Mode**: Advanced analysis for comparing haplotype profiles
- **Docker Support**: Easy deployment with pre-configured Docker images

## Installation

### Docker images (preferred)
The preferred way to run MatchY is using the provided Docker images. This ensures that all dependencies are correctly installed and configured.

#### GUI (Streamlit Dashboard)
1. Install Docker (Desktop) from https://www.docker.com/get-started
2. Pull the Docker image: `docker pull dionzand/matchy:latest`
3. Run the Docker container: `docker run -p 8501:8501 dionzand/matchy:latest`
4. Access the Streamlit dashboard at `http://localhost:8501`

#### CLI (Command Line Interface)
1. Install Docker (Desktop) from https://www.docker.com/get-started
2. Pull the Docker image: `docker pull dionzand/match-cli:latest`
3. Run the CLI:
   ```bash
   # Show help
   docker run --rm dionzand/match-cli:latest

   # Run with a config file (mount your config directory)
   docker run --rm -v /path/to/your/data:/app/data dionzand/match-cli:latest --config data/config.ini

   # Run with specific options
   docker run --rm -v /path/to/your/data:/app/data dionzand/match-cli:latest --config data/config.ini --skip-inside

   # Run trace mode
   docker run --rm -v /path/to/your/data:/app/data dionzand/match-cli:latest --config data/config.ini --trace-mode
   ```

### Local installation (Linux only)
If you prefer to run MatchY locally, follow these steps. **Note: Local installation is only supported on Linux.**

1. Install python
2. Clone the repository: `git clone https://github.com/genid/MatchY.git`
3. Install dependencies: `pip install -r requirements.txt`
4. Run the application
   * CLI interface: `python main.py`
   * Streamlit Dashboard: `streamlit run streamlit_app.py`

## Quick Start

### Using the GUI

1. Launch the Streamlit dashboard (see [Installation](#installation))
2. Navigate through the pages:
   - **Home**: Overview and software settings
   - **Pedigree Builder**: Create and visualize family trees
   - **Haplotype Editor**: Input Y-STR haplotype data
   - **Simulation**: Configure and run simulations
   - **Results**: View and export analysis results

### Using the CLI

```bash
# Basic simulation
python main.py --config config.ini

# Skip inside/outside pedigree calculations
python main.py --config config.ini --skip-inside
python main.py --config config.ini --skip-outside

# Run in trace mode
python main.py --config config.ini --trace-mode
```

## Documentation

- **[User Manual](USER_MANUAL.md)**: Comprehensive guide for using MatchY
- **[Parameters Reference](PARAMETERS_REFERENCE.md)**: Detailed description of all configuration parameters
- **[Contributing Guide](CONTRIBUTING.md)**: Guidelines for contributors

## Contributing

We welcome contributions from the community! Please see our [Contributing Guide](CONTRIBUTING.md) for details on:

- Setting up your development environment
- Coding standards and style guidelines
- Submitting pull requests
- Reporting bugs and suggesting features

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2026 Department of Pathology and Clinical Bioinformatics, Erasmus MC University Medical Center Rotterdam, The Netherlands; Institute of Medical Informatics and Statistics, Kiel University, University Hospital Schleswig-Holstein, Kiel, Germany; Chair of Epidemiology, Medical Biometry and Medical Informatics, Department of Medicine, Health and Medical University Erfurt, Erfurt, Germany.

## Citation

If you use MatchY in your research, please cite:

```
[Citation information to be added]
```

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/genid/MatchY/issues)
- **Questions**: Check existing issues or create a new one with the "question" label
- **Documentation**: Refer to the [User Manual](USER_MANUAL.md) for detailed usage instructions

---

**Developed for the forensic genetics community**