# pedigreeLR
pedigreeLR is a powerful pedigree-based tool designed to estimate match probabilities for Y-chromosomal haplotypes. Its mathematical framework leverages marker mutation rates, pedigree structure, and the known haplotypes of individuals within the family tree. By combining this data, the tool accurately estimates match probabilities with a person of interest using a Monte Carlo simulation to model mutations.

To enhance efficiency, pedigreeLR employs importance sampling, significantly streamlining the simulation process. It supports any number of Y-STR markers, including multi-copy markers and intermediate alleles.

## Installation

1. Install python
2. Install dependencies: `pip install -r requirements.txt`
3. Run the application
   * CLI interface: `python main.py`
   * Streamlit Dashboard: `streamlit run streamlit_app.py`

## Formatting

To format the code automatically, use the following commands:

```
ruff format pedigree
ruff check --fix --select I
```

Make sure to have ruff installed (`pip install ruff`)

## Type checking

To type check the code for inconsistencies you can use mypy:

```
mypy pedigree
```

Make sure to have mypy installed (`pip install mypy`)
