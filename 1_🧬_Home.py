from __future__ import annotations
from configparser import ConfigParser
from datetime import datetime
from pathlib import Path
from io import StringIO
import streamlit as st
import pandas as pd
from pedigree_lr.data import load_pedigree_from_upload, get_marker_set_names, load_marker_set_from_database
from pedigree_lr.models import SimulationResult, SimulationParameters
from pedigree_lr.reporting import StreamlitReporter, create_html_pdf_report, setup_logger_streamlit
from pedigree_lr.simulation import run_simulation
from pedigree_lr.visualization import st_visualize_pedigree

st.set_page_config(
    page_title="MatchY",
    page_icon="icon.png",
    layout="wide",
    initial_sidebar_state="expanded",
)

logger = setup_logger_streamlit()

st.logo("logo_minimal.png", icon_image="icon.png")
st.markdown(body=
            '''
            <style>
            /* Default size when sidebar is open */
                section[data-testid="stSidebar"][aria-expanded="true"] img[data-testid="stSidebarLogo"] {
                  height: 70px; /* or whatever height you want */
                  margin-top: 0.75rem;
                  transition: height 0.3s ease;
                }
                
                /* Smaller size when sidebar is closed */
                section[data-testid="stSidebar"][aria-expanded="false"] img[data-testid="stLogo"] {
                  height: 50px; /* smaller logo */
                  transition: height 0.3s ease;
                }
            </style>
            ''', unsafe_allow_html=True)

config_path = Path(__file__).resolve().parent / "data" / "config.ini"
global_config = ConfigParser()
global_config.optionxform = str  # type: ignore
global_config.read(config_path)
st.session_state.global_config = global_config


def render_simulation() -> SimulationResult | None:
    st_visualize_pedigree(st.session_state.pedigree,
                          global_config=st.session_state.global_config)

    input_placeholder = st.empty()
    with input_placeholder.container():
        possible_suspects = [individual.name for individual in st.session_state.pedigree.individuals
                             if individual.haplotype_class == "known" or individual.haplotype_class == "suspect"]

        if len(possible_suspects) == 0:
            st.warning(
                "Make sure the pedigree contains at least one individual with known or haplotype to be able to select a suspect.")
            return None

        col1, col2 = st.columns(2)

        suspect_name = col1.selectbox(
            "Select suspect/trace",
            options=possible_suspects,
            index=0,
            help="Select which known individual is the suspect/trace.",
        )

        possible_excluded_individuals = [individual.name for individual in st.session_state.pedigree.individuals
                                         if individual.name != suspect_name]

        excluded_individuals = col2.multiselect(
            "Choose individuals to exclude from the simulation",
            options=possible_excluded_individuals,
            help="Exclude individuals based on case context information.",
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

            two_step_mutation_factor = col3.number_input(
                "Two-step mutation factor",
                min_value=0.0,
                max_value=1.0,
                value=global_config.getfloat("simulation_parameters", "two_step_mutation_factor"),
                step=0.01,
                format="%.2f",
                help="The factor by which the mutation rate is multiplied for two-step mutations. "
                     "This is used to simulate the effect of two-step mutations on the probability of a match.",
            )

            stability_window = col3.number_input(
                "Stability window",
                min_value=0,
                value=global_config.getint("simulation_parameters", "stability_window"),
                step=1,
                help="The number of iterations the simulation should be stable.",
            )

            model_validity = col3.number_input(
                "Model validity",
                min_value=0.0,
                value=global_config.getfloat("simulation_parameters", "model_validity"),
                step=0.0001,
                format="%.4f",
                help="The maximum allowed relative difference between results for the model to be considered valid.",
            )

            skip_inside = col4.checkbox(
                "Skip inside pedigree probabilities",
                value=False,
                help="If checked, the simulation will skip calculating inside pedigree probabilities.",
            )

            skip_outside = col4.checkbox(
                "Skip outside pedigree probabilities",
                value=False,
                help="If checked, the simulation will skip calculating outside pedigree probabilities.",
            )

            trace_mode = col4.checkbox(
                "Use trace mode",
                value=False,
                help="If checked, expects a trace instead of a suspect.",
            )

        simulation_name = st.text_input("Give this simulation a name").replace(" ", "_").lower()
        user_name = st.text_input("Your name", help="Enter your name to be included in the report.")

        if not st.button("Start simulation",
                         type="primary", ):
            return None

    input_placeholder.empty()

    if st.button("Stop simulation and reset parameters",
                 type="secondary", ):
        st.rerun()

    st.warning("Simulation in progress. Do not close the browser or navigate away from this page...")

    progress_placeholder = st.empty()
    log_placeholder = st.empty()

    reporter = StreamlitReporter(
        log_container=log_placeholder.container(height=600, border=True),
        progress_container=progress_placeholder.container()
    )

    folder_name = f"{simulation_name.replace(' ', '_').lower()}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    results_path = Path(__file__).resolve().parent / "results" / folder_name
    results_path.mkdir(parents=True, exist_ok=True)

    simulation_parameters = SimulationParameters(
        two_step_mutation_factor=two_step_mutation_factor,
        stability_window=stability_window,
        model_validity_threshold=model_validity,
        simulation_name=simulation_name,
        number_of_threads=1,
        results_path=results_path,
        user_name=user_name,
    )

    simulation_result = run_simulation(
        input_pedigree=st.session_state.pedigree,
        suspect_name=st.session_state.suspect,
        marker_set=st.session_state.marker_set,
        simulation_parameters=simulation_parameters,
        reporter=reporter,
        skip_inside=skip_inside,
        skip_outside=skip_outside,
        trace_mode=trace_mode
    )

    progress_placeholder.empty()

    st.text("Results:")

    try:
        proposal_distribution_dataframe = pd.DataFrame(list(simulation_result.inside_match_probability.items()),
                                                       columns=['Number of matches', 'Probability'])
        proposal_distribution_dataframe.set_index('Number of matches', inplace=True)
        st.dataframe(proposal_distribution_dataframe.sort_values(by='Number of matches'))
    except Exception as e:
        st.error(f"Error displaying proposal distribution: {e}")

    st.text("Outside match probability:")
    st.write(f"{simulation_result.outside_match_probability:.4f}")

    report_bytes = create_html_pdf_report(simulation_result)
    st.download_button(
        label="Download simulation report",
        data=report_bytes,
        file_name=f"{simulation_result.simulation_parameters.simulation_name}_report.pdf",
        mime="application/pdf",
        help="Download the simulation report as a PDF file."
    )

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
                                         help="Upload a pedigree file in TGF or PED format. Make sure the node labels correspond to the names in the haplotypes file.",
                                         accept_multiple_files=False)

        haplotypes_file = st.file_uploader("Upload haplotypes file",
                                           type=["json"],
                                           help="Upload haplotypes file in JSON format. Use Haplotype editor to create a haplotype file.",
                                           accept_multiple_files=False)

        if st.button("Upload files",
                     type="primary"):
            if selected_marker_set:
                st.session_state.marker_set = load_marker_set_from_database(
                    selected_marker_set)
            if pedigree_file is not None:
                file_extension = Path(pedigree_file.name).suffix
                stringio = StringIO(pedigree_file.getvalue().decode("utf-8"))
                st.session_state.pedigree = load_pedigree_from_upload(stringio, file_extension)
                if st.session_state.pedigree is not None:
                    st.success("Pedigree file uploaded successfully")
            if haplotypes_file is not None:
                stringio = StringIO(haplotypes_file.getvalue().decode("utf-8"))
                st.session_state.pedigree.read_known_haplotypes_from_file(stringio, st.session_state.marker_set)
                st.success("Haplotypes file uploaded successfully")

        st.divider()

        with open(r"examples/example.tgf") as file:
            st.download_button("Download example pedigree file in TGF format", file, "example_pedigree.tgf",
                               type="tertiary")

    if st.session_state.pedigree is not None:
        result = render_simulation()
    else:
        st.error(
            "Please upload all necessary files via the sidebar, or enter them via Pedigree builder and Haplotypes editor.")
