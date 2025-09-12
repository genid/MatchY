# MatchY
MatchY is a powerful pedigree-based tool designed to estimate match probabilities for Y-STR haplotypes. Its mathematical framework leverages marker mutation rates, pedigree structure, and the known haplotypes of individuals within the family tree. By combining this data, the tool accurately estimates match probabilities with a person of interest using a Monte Carlo simulation with importance sampling to model mutations.

MatchY supports any number of Y-STR markers, including multi-copy markers and intermediate alleles.

## Installation

### Docker image (preferred)
The preferred way to run MatchY is using the provided Docker image. This ensures that all dependencies are correctly installed and configured.

1. Install Docker (Desktop) from https://www.docker.com/get-started
2. Pull the Docker image: `docker pull dionzand/matchy:latest`
3. Run the Docker container: `docker run -p 8501:8501 dionzand/matchy:latest`
4. Access the Streamlit dashboard at `http://localhost:8501`

### Local installation
If you prefer to run MatchY locally, follow these steps:

1. Install python
2. Clone the repository: `git clone https://github.com/genid/MatchY.git`
3. Install dependencies: `pip install -r requirements.txt`
3. Run the application
   * CLI interface: `python main.py`
   * Streamlit Dashboard: `streamlit run streamlit_app.py`