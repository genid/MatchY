from pathlib import Path
from random import Random

import streamlit as st

from pedigree.config import load_config
from pedigree.data import load_marker_set, load_pedigree
from pedigree.simulation import run_simulation
from pedigree.visualization import st_visualize_pedigree, st_visualize_simulation_result

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

    graph_placeholder = st.empty()
    with graph_placeholder:
        st_visualize_pedigree(pedigree)

    progress_bar = st.progress(0, text="Simulating")

    def _update_progress(count: int, total: int) -> None:
        progress_bar.progress(count / total, text="Simulating")

    result = run_simulation(
        pedigree=pedigree,
        marker_set=marker_set,
        suspect=config.suspect,
        number_of_iterations=config.number_of_iterations,
        random=Random(config.random_seed),
        show_progress=lambda count, total: progress_bar.progress(count / total, text="Simulating")
    )

    progress_bar.empty()

    st_visualize_simulation_result(result, config.number_of_iterations)
