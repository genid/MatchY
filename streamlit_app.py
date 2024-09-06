from pathlib import Path
from random import Random

import pandas as pd
import streamlit as st
from streamlit.delta_generator import DeltaGenerator

from pedigree_lr.config import load_config
from pedigree_lr.data import load_marker_set, load_pedigree
from pedigree_lr.models import SimulationResult
from pedigree_lr.reporting import StreamlitReporter
from pedigree_lr.simulation import run_simulation
from pedigree_lr.visualization import st_visualize_pedigree

_data_dir = Path("data")

st.set_page_config(
        page_title="fraternitY",
        page_icon="🧬",
        layout="wide"
    )


def render_simulation() -> SimulationResult | None:
    input_placeholder = st.empty()

    with input_placeholder.container():
        col1, col2 = st.columns(2)

        number_of_iterations = col1.number_input(
            "Number of iterations",
            min_value=1000,
            max_value=10000000,
            value=config.number_of_iterations,
            step=1000,
        )

        random_seed = col2.number_input(
            "Random seed",
            min_value=0,
            max_value=1000000,
            value=config.random_seed,
            step=1,
        )

        if not st.button("Start simulation"):
            return None

    input_placeholder.empty()

    progress_placeholder = st.empty()
    log_placeholder = st.empty()

    reporter = StreamlitReporter(
        log_container=log_placeholder.container(height=600, border=True),
        progress_container=progress_placeholder.container()
    )

    result = run_simulation(
        pedigree=pedigree,
        suspect_name=config.suspect,
        marker_set=marker_set,
        number_of_iterations=number_of_iterations,
        random=Random(random_seed),
        reporter=reporter,
        show_simulated_pedigrees=config.show_simulated_pedigrees
    )

    progress_placeholder.empty()
    log_placeholder.empty()

    st.text(f"run_time_pedigree_probability: {result.run_time_pedigree_probability}")
    st.text(f"run_time_proposal_distribution: {result.run_time_proposal_distribution}")

    st.text(f"average_pedigree_probability: {result.average_pedigree_probability}")

    st.text(f"proposal_distribution:")
    st.table(
        pd.DataFrame(
            result.proposal_distribution.values(),
            index=list(result.proposal_distribution.keys()),
            columns=["probability"]
        )
    )

    if st.button("New simulation"):
        st.rerun()

    return result


if __name__ == '__main__':
    config = load_config(Path("config.ini"))
    marker_set = load_marker_set(config)
    pedigree = load_pedigree(config, marker_set)

    col1, col2 = st.columns(2)

    graph_placeholder = col1.empty()
    with graph_placeholder:
        st_visualize_pedigree(pedigree)

    with col2:
        result = render_simulation()
