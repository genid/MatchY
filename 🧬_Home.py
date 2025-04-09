from __future__ import annotations
from datetime import datetime
from pathlib import Path
from random import Random
from io import StringIO
import streamlit as st
import pandas as pd
from pedigree_lr.data import load_pedigree_from_upload, get_marker_set_names, load_marker_set_from_database
from pedigree_lr.models import SimulationResult
from pedigree_lr.reporting import StreamlitReporter
from pedigree_lr.simulation import run_simulation
from pedigree_lr.visualization import st_visualize_pedigree

st.set_page_config(
    page_title="match-Y",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def render_simulation() -> SimulationResult | None:
    input_placeholder = st.empty()

    with input_placeholder.container():
        st_visualize_pedigree(st.session_state.pedigree)

        possible_suspects = [individual.name for individual in st.session_state.pedigree.individuals
                             if individual.haplotype_class == "known" or individual.haplotype_class == "suspect"]

        col1, col2 = st.columns(2)

        suspect_name = col1.selectbox(
            "Select suspect",
            options=possible_suspects,
            index=0,
            help="Select the individual you want to test against the pedigree.",
        )

        possible_excluded_individuals = [individual.name for individual in st.session_state.pedigree.individuals
                                         if individual.name != suspect_name]

        excluded_individuals = col2.multiselect(
            "Choose individuals to exclude from the simulation",
            options=possible_excluded_individuals,
            help="Choose individuals to exclude from the simulation. "
                 "The excluded individuals must be in the pedigree.",
        )

        if st.button("Set suspect and unknown individuals"):
            st.session_state.suspect = suspect_name
            st.session_state.pedigree.set_suspect(st.session_state.suspect)
            st.session_state.pedigree.exclude_individuals(excluded_individuals)
            st.rerun()

        if st.session_state.get("suspect", None) is None:
            return None

        with st.expander("Set simulation parameters", expanded=True):
            col3, col4 = st.columns(2)

            number_of_iterations = col3.number_input(
                "Number of iterations",
                min_value=1000,
                value=1000000,
                step=1000,
            )

            two_step_mutation_factor = col3.number_input(
                "Two-step mutation factor",
                min_value=0.0,
                max_value=1.0,
                value=0.03,
                step=0.01,
                format="%.2f",
                help="The factor by which the mutation rate is multiplied for two-step mutations. "
                     "This is used to simulate the effect of two-step mutations on the probability of a match.",
            )

            stability_window = col3.number_input(
                "Stability window",
                min_value=0,
                value=1000,
                step=1,
                help="The number of iterations the simulation should be stable.",
            )

            stability_min_iterations = col3.number_input(
                "Stability minimum iterations",
                min_value=0,
                value=2000,
                step=1,
                help="The minimum number of iterations the simulation should run before testing for stability.",
            )

            random_seed = col4.number_input(
                "Random seed",
                min_value=0,
                value=1234,
                step=1,
            )

            stability_threshold = col4.number_input(
                "Stability threshold",
                min_value=0.0000,
                value=0.0001,
                step=0.0001,
                format="%.4f",
                help="The maximum allowed relative change between consecutive probabilities for stability.",
            )

            model_validity = col4.number_input(
                "Model validity",
                min_value=0.0,
                value=0.05,
                step=0.0001,
                format="%.4f",
                help="The maximum allowed relative difference between results for the model to be considered valid.",
            )

        if not st.button("Start simulation",
                         type="primary", ):
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
    if st.download_button("Download results as report",
                          simulation_result.download_results(random_seed=random_seed),
                          f"matchY_report_{datetime.now().strftime('%Y%m%d%H%M%S')}.txt",
                          type="primary", ):
        st.success("Download started")

    try:
        proposal_distribution_dataframe = pd.DataFrame(list(simulation_result.proposal_distribution.items()),
                                                       columns=['Number of matches', 'Probability'])
        proposal_distribution_dataframe.set_index('Number of matches', inplace=True)
        st.dataframe(proposal_distribution_dataframe.sort_values(by='Number of matches'))
    except Exception as e:
        st.error(f"Error displaying proposal distribution: {e}")

    st.text("Outside match probability:")
    st.write(f"{simulation_result.outside_match_probability:.4f}")

    if st.button("New simulation"):
        st.rerun()

    return simulation_result


if __name__ == '__main__':
    if "marker_set" not in st.session_state:
        st.session_state.marker_set = None
    if "pedigree" not in st.session_state:
        st.session_state.pedigree = None

    with st.sidebar:
        selected_marker_set = st.selectbox("Select marker set",
                                           get_marker_set_names(),
                                           index=0,
                                           help="Select the marker set you want to use for the simulation. "
                                                "The marker set defines the markers and their mutation rates.",
                                           )

        pedigree_file = st.file_uploader("Upload pedigree file",
                                         type=["tgf", "ped"],
                                         help="Upload a pedigree file in TGF or PED format. Make sure the node labels correspond to the haplotype file names.",
                                         accept_multiple_files=False)

        haplotypes_files = st.file_uploader("Upload haplotypes file(s)",
                                            type=["csv"],
                                            help="Upload haplotypes file(s). File names are used as individual names.",
                                            accept_multiple_files=True)

        st.divider()

        if st.button("Upload files",
                     type="primary", ):
            if selected_marker_set:
                st.session_state.marker_set = load_marker_set_from_database(
                    selected_marker_set)
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

            st.success("Files uploaded successfully")

        st.divider()

        with open(r"examples/RM/mutation_rates.csv") as file:
            st.download_button("Download example marker set file", file, "mutation_rates.csv")
        with open(r"examples/pedigree_large.tgf") as file:
            st.download_button("Download example pedigree file in TGF format", file, "pedigree_large.tgf")
        with open(r"examples/RM/George.csv") as file:
            st.download_button("Download example haplotypes file", file, "George.csv")

    if st.session_state.marker_set is not None and st.session_state.pedigree is not None:
        result = render_simulation()
    else:
        st.error("Please upload all necessary files via the sidebar")
