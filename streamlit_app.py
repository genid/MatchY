from __future__ import annotations

from pathlib import Path
from random import Random

from io import StringIO

import pandas as pd
import streamlit as st

from pedigree_lr.data import load_marker_set_from_upload, load_pedigree_from_upload
from pedigree_lr.models import SimulationResult
from pedigree_lr.reporting import StreamlitReporter
from pedigree_lr.simulation import run_simulation
from pedigree_lr.visualization import st_visualize_pedigree

_data_dir = Path("data")

st.set_page_config(
        page_title="pedigreeLR",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def render_simulation() -> SimulationResult | None:
    input_placeholder = st.empty()

    with input_placeholder.container():
        col1, col2 = st.columns(2)

        number_of_iterations = col1.number_input(
            "Number of iterations",
            min_value=1000,
            max_value=10000000,
            value=10000,
            step=1000,
        )

        random_seed = col2.number_input(
            "Random seed",
            min_value=0,
            max_value=1000000,
            value=1234,
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
        pedigree=st.session_state.pedigree,
        suspect_name=st.session_state.suspect,
        marker_set=st.session_state.marker_set,
        number_of_iterations=number_of_iterations,
        random=Random(random_seed),
        reporter=reporter,
        show_simulated_pedigrees=False
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
    if "marker_set" not in st.session_state:
        st.session_state.marker_set = None
        st.error("Please upload all necessary files via the sidebar")
        st.markdown("Upload a marker set file, a pedigree file, and haplotype file(s) via the sidebar.\n"
                    "You can alter the pedigree file in TGF format using the yEd software.\n"
                    "Make sure the marker set file and haplotype files all contain headers.\n"
                    "All markers in the marker set file should also be present in the haplotype files.\n"
                    "Make sure that all marker and individuals names are identical between files.\n")

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

        st.session_state.suspect = st.text_input("Suspect")
        st.warning("The suspect must be in the pedigree")

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
                    st.session_state.pedigree.reroot_pedigree(st.session_state.suspect)

    if st.session_state.pedigree is not None:
        st_visualize_pedigree(st.session_state.pedigree)

    if st.session_state.suspect != "":
        result = render_simulation()
