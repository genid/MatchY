# MatchY
MatchY is a powerful pedigree-based tool designed to estimate match probabilities for Y-STR haplotypes. Its mathematical framework leverages marker mutation rates, pedigree structure, and the known haplotypes of individuals within the family tree. By combining this data, the tool accurately estimates match probabilities with a person of interest using a Monte Carlo simulation with importance sampling to model mutations.

MatchY supports any number of Y-STR markers, including multi-copy markers and intermediate alleles.

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