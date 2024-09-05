from pathlib import Path
from random import Random

import pandas as pd
import streamlit as st

from pedigree.config import load_config
from pedigree.data import load_marker_set, load_pedigree
from pedigree.reporting import StreamlitReporter
from pedigree.simulation import run_simulation
from pedigree.visualization import st_visualize_pedigree

_data_dir = Path("data")

st.set_page_config(
        page_title="fraternitY",
        page_icon="🧬",
        layout="wide"
    )


if __name__ == '__main__':
    config = load_config(Path("config.ini"))
    marker_set = load_marker_set(config)
    pedigree = load_pedigree(config, marker_set)

    col1, col2 = st.columns(2)

    graph_placeholder = col1.empty()
    with graph_placeholder:
        st_visualize_pedigree(pedigree)

    progress_placeholder = col2.empty()
    log_placeholder = col2.empty()

    reporter = StreamlitReporter(
        log_container=log_placeholder.container(height=600, border=True),
        progress_container=progress_placeholder.container()
    )

    result = run_simulation(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect_name=config.suspect,
        number_of_iterations=config.number_of_iterations,
        random=Random(config.random_seed),
        reporter=reporter,
    )

    progress_placeholder.empty()
    log_placeholder.empty()

    col2.text(f"run_time_pedigree_probability: {result.run_time_pedigree_probability}")
    col2.text(f"run_time_proposal_distribution: {result.run_time_proposal_distribution}")

    col2.text(f"average_pedigree_probability: {result.average_pedigree_probability}")

    col2.text(f"proposal_distribution:")
    col2.table(
        pd.DataFrame(
            result.proposal_distribution.values(),
            index=list(result.proposal_distribution.keys()),
            columns=["probability"]
        )
    )
