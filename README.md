# Pedigree LR

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
