from __future__ import annotations

from inspect import stack
from pathlib import Path
from random import Random
from io import StringIO
import streamlit as st
import pandas as pd
from pedigree_lr.data import load_marker_set_from_upload, load_pedigree_from_upload
from pedigree_lr.models import SimulationResult
from pedigree_lr.reporting import StreamlitReporter
from pedigree_lr.simulation import run_simulation
from pedigree_lr.visualization import st_visualize_pedigree

_data_dir = Path("data")

st.set_page_config(
        page_title="match-Y",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def render_simulation() -> SimulationResult | None:
    input_placeholder = st.empty()

    with input_placeholder.container():
        with st.expander("Set simulation parameters", expanded=False):
            col1, col2 = st.columns(2)

            number_of_iterations = col1.number_input(
                "Number of iterations",
                min_value=1000,
                max_value=10000000,
                value=10000,
                step=1000,
            )

            two_step_mutation_factor = col1.number_input(
                "Two-step mutation factor",
                min_value=0.0,
                max_value=1.0,
                value=0.03,
                step=0.01,
                format="%.2f",
                help="The factor by which the mutation rate is multiplied for two-step mutations. "
                     "This is used to simulate the effect of two-step mutations on the probability of a match.",
            )

            stability_window = col1.number_input(
                "Stability window",
                min_value=0,
                value=500,
                step=1,
                help="The number of iterations the simulation should be stable.",
            )

            stability_min_iterations = col1.number_input(
                "Stability minimum iterations",
                min_value=0,
                value=2000,
                step=1,
                help="The minimum number of iterations the simulation should run before testing for stability.",
            )

            random_seed = col2.number_input(
                "Random seed",
                min_value=0,
                max_value=1000000,
                value=1234,
                step=1,
            )

            stability_threshold = col2.number_input(
                "Stability threshold",
                min_value=0.0000,
                value=0.0001,
                step=0.0001,
                format="%.4f",
                help="The maximum allowed relative change between consecutive probabilities for stability.",
            )

            model_validity = col2.number_input(
                "Model validity",
                min_value=0.0,
                value=0.005,
                step=0.0001,
                format="%.4f",
                help="The maximum allowed relative difference between results for the model to be considered valid.",
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

    simulation_parameters = {
        "number_of_iterations": number_of_iterations,
        "two_step_mutation_factor": two_step_mutation_factor,
        "stability_window": stability_window,
        "stability_min_iterations": stability_min_iterations,
        "stability_threshold": stability_threshold,
        "model_validity_threshold": model_validity
    }

    simulation_result = run_simulation(
        pedigree=st.session_state.pedigree,
        suspect_name=st.session_state.suspect,
        marker_set=st.session_state.marker_set,
        simulation_parameters=simulation_parameters,
        random=Random(random_seed),
        reporter=reporter,
    )

    progress_placeholder.empty()

    st.text("Results:")
    proposal_distribution_dataframe = pd.DataFrame(list(simulation_result.proposal_distribution.items()), columns=['Number of matches', 'Probability'])
    proposal_distribution_dataframe.set_index('Number of matches', inplace=True)
    st.dataframe(proposal_distribution_dataframe.sort_values(by='Number of matches'))

    st.text("Outside match probability:")
    st.write(f"{simulation_result.outside_match_probability:.4f}")

    if st.download_button("Download results as report", simulation_result.download_results(random_seed=random_seed), "matchY_report.txt"):
        st.success("Download started")

    if st.button("New simulation"):
        st.rerun()

    return simulation_result


if __name__ == '__main__':
    if "marker_set" not in st.session_state:
        st.session_state.marker_set = None
        st.error("Please upload all necessary files via the sidebar")

    if "pedigree" not in st.session_state:
        st.session_state.pedigree = None
    if "suspect" not in st.session_state:
        st.session_state.suspect = None

    with st.sidebar:
        marker_set_file = st.file_uploader("Upload marker set file",
                                           type=["csv"],
                                           help="Upload a marker set file in CSV format, with columns 'name' and 'mutation_rate'. File should contain a header.",
                                           accept_multiple_files=False)

        with open(r"examples/RM/mutation_rates.csv") as file:
            st.download_button("Download example marker set file", file, "mutation_rates.csv")

        pedigree_file = st.file_uploader("Upload pedigree file",
                                         type=["tgf", "ped"],
                                         help="Upload a pedigree file in TGF or PED format. Make sure the node labels correspond to the haplotype file names.",
                                         accept_multiple_files=False)

        with open(r"examples/pedigree_large.tgf") as file:
            st.download_button("Download example pedigree file in TGF format", file, "pedigree_large.tgf")

        haplotypes_files = st.file_uploader("Upload haplotypes file(s)",
                                            type=["csv"],
                                            help="Upload haplotypes file(s). File names are used as individual names.",
                                            accept_multiple_files=True)

        with open(r"examples/RM/George.csv") as file:
            st.download_button("Download example haplotypes file", file, "George.csv")

        st.session_state.suspect = st.text_input("Enter suspect name")
        st.warning("The suspect must be in the pedigree")

        excluded_individuals = st.text_input("Enter excluded individuals (comma separated)")
        st.session_state.excluded_individuals = [ind.strip() for ind in excluded_individuals.split(",")]
        st.warning("The excluded individuals must be in the pedigree")

        upload_files = st.button("Upload files")
        if upload_files:
            if marker_set_file is not None:
                st.session_state.marker_set = load_marker_set_from_upload(
                    StringIO(marker_set_file.getvalue().decode("utf-8")))
            if pedigree_file is not None:
                file_extension = Path(pedigree_file.name).suffix
                stringio = StringIO(pedigree_file.getvalue().decode("utf-8"))
                st.session_state.pedigree = load_pedigree_from_upload(stringio, file_extension)

            if haplotypes_files is not None:
                for haplotypes_file in haplotypes_files:
                    stringio = StringIO(haplotypes_file.getvalue().decode("utf-8"))
                    name = haplotypes_file.name.split(".")[0]
                    st.session_state.pedigree.read_known_haplotype_from_file(name, stringio,
                                                                             st.session_state.marker_set)

            st.session_state.pedigree.exclude_individuals(st.session_state.excluded_individuals)

    if st.session_state.pedigree is not None:
        st_visualize_pedigree(st.session_state.pedigree)

    if st.session_state.suspect != "":
        result = render_simulation()
