from __future__ import annotations
from configparser import ConfigParser
from datetime import datetime
from pathlib import Path
from io import StringIO
import multiprocessing
import streamlit as st
import pandas as pd
from pedigree_lr.data import load_pedigree_from_upload, get_marker_set_names, load_marker_set_from_database
from pedigree_lr.models import SimulationResult, SimulationParameters
from pedigree_lr.reporting import StreamlitReporter, create_html_pdf_report, create_trace_mode_report, setup_logger_streamlit
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
                "⚠️ The pedigree must contain at least one individual with known haplotype to select as a suspect.")
            return None

        st.markdown("### ⚙️ Simulation Configuration")

        col1, col2 = st.columns(2)

        # Check if TRACE profile is available
        has_trace_profile = st.session_state.get("trace", None) is not None

        # Trace mode checkbox (only show if TRACE profile is available)
        trace_mode = False
        if has_trace_profile:
            # Use loaded config if available
            default_trace_mode = False
            if "loaded_config" in st.session_state and st.session_state.loaded_config:
                default_trace_mode = st.session_state.loaded_config.get("trace_mode", False)

            trace_mode = col1.checkbox(
                "Use trace mode",
                value=default_trace_mode,
                help="If checked, uses the TRACE profile from the haplotypes JSON file instead of a suspect.",
            )

            # Trace status validation
            if trace_mode:
                st.success("✅ TRACE profile loaded and active")
            else:
                st.info("ℹ️ TRACE profile detected. Enable 'Use trace mode' to activate it.")

        # Suspect selection (disabled when trace mode active)
        # Use loaded config for default suspect
        default_suspect_index = 0
        if "loaded_config" in st.session_state and st.session_state.loaded_config:
            loaded_suspect = st.session_state.loaded_config.get("suspect", None)
            if loaded_suspect and loaded_suspect in possible_suspects:
                default_suspect_index = possible_suspects.index(loaded_suspect)

        suspect_name = col1.selectbox(
            "Select suspect",
            options=possible_suspects,
            index=default_suspect_index,
            help="Select which known individual is the suspect." if not trace_mode else "Disabled in trace mode - using TRACE profile.",
            disabled=trace_mode,
        )

        possible_excluded_individuals = [individual.name for individual in st.session_state.pedigree.individuals
                                         if individual.name != suspect_name]

        # Use loaded config for default excluded individuals
        default_excluded = []
        if "loaded_config" in st.session_state and st.session_state.loaded_config:
            loaded_excluded = st.session_state.loaded_config.get("exclude_individuals", [])
            # Filter to only include individuals that exist in possible_excluded_individuals
            default_excluded = [ind for ind in loaded_excluded if ind in possible_excluded_individuals and ind]

        excluded_individuals = col2.multiselect(
            "Choose individuals to exclude from the simulation",
            options=possible_excluded_individuals,
            default=default_excluded,
            help="Exclude individuals based on case context information.",
        )

        if st.button("✅ Confirm Selection", type="primary"):
            if not trace_mode:
                st.session_state.suspect = suspect_name
                st.session_state.pedigree.set_suspect(st.session_state.suspect)
            st.session_state.pedigree.exclude_individuals(excluded_individuals)
            st.rerun()

        if not trace_mode and st.session_state.get("suspect", None) is None:
            return None

        st.markdown("")

        with st.expander("⚙️ Simulation Parameters", expanded=True):
            col3, col4 = st.columns(2)

            # Use loaded config values if available
            default_two_step = global_config.getfloat("simulation_parameters", "two_step_mutation_fraction")
            default_stability = global_config.getint("simulation_parameters", "batch_length")
            default_validity = global_config.getfloat("simulation_parameters", "convergence_criterion")
            default_adaptive_bias = False

            if "loaded_config" in st.session_state and st.session_state.loaded_config:
                default_two_step = st.session_state.loaded_config.get("two_step_mutation_fraction", default_two_step)
                default_stability = st.session_state.loaded_config.get("batch_length", default_stability)
                default_validity = st.session_state.loaded_config.get("convergence_criterion", default_validity)
                # Note: We can't determine from config alone if adaptive bias was used
                # (both adaptive bias and default bias result in bias=None in config)
                # So we default to False and let user check it manually if needed

            two_step_mutation_factor = col3.number_input(
                "Two-step mutation fraction",
                min_value=0.0,
                max_value=1.0,
                value=default_two_step,
                step=0.01,
                format="%.2f",
                help="The fraction by which the mutation rate is multiplied for two-step mutations. "
                     "This is used to simulate the effect of two-step mutations on the probability of a match.",
            )

            stability_window = col3.number_input(
                "Batch length",
                min_value=0,
                value=default_stability,
                step=1,
                help="The number of iterations the simulation should be stable.",
            )

            model_validity = col3.number_input(
                "Convergence criterion",
                min_value=0.0,
                value=default_validity,
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

            adaptive_bias = col4.checkbox(
                "Use adaptive bias mode",
                value=default_adaptive_bias,
                help="If checked, dynamically adjusts bias values for each model based on performance. "
                     "Best-performing models get lower bias (0.05), worst get higher bias (0.25).",
            )

            # Bias value input (only enabled when adaptive bias is not selected)
            # Get default bias value from loaded config or None
            default_bias = None
            if "loaded_config" in st.session_state and st.session_state.loaded_config:
                default_bias = st.session_state.loaded_config.get("bias", None)

            # Override to None if adaptive bias is enabled
            if adaptive_bias:
                bias_value = None
            else:
                bias_value = col4.number_input(
                    "Bias value",
                    min_value=0.0,
                    max_value=1.0,
                    value=default_bias,
                    step=0.01,
                    format="%.2f",
                    help="Explicit bias value for importance sampling. Leave empty for default bias mode (dynamically calculated based on pedigree structure). Only used when adaptive bias mode is disabled.",
                    placeholder="Default (auto)",
                )

        st.markdown("")

        # Execution mode selection
        with st.expander("🚀 Execution Mode", expanded=True):
            # Determine default number of threads from loaded config
            default_n_threads = min(4, multiprocessing.cpu_count())
            if "loaded_config" in st.session_state and st.session_state.loaded_config:
                loaded_threads = st.session_state.loaded_config.get("number_of_threads", default_n_threads)
                default_n_threads = min(loaded_threads, multiprocessing.cpu_count())

            execution_mode = st.radio(
                "Choose how to run the simulation:",
                options=["Run in Browser (Single-threaded)", "Submit to CLI (Multi-threaded)"],
                help=(
                    "Run in Browser: Simpler, but slower (single-threaded). "
                    "Submit to CLI: Faster with multiprocessing, but creates temporary files."
                )
            )

            if execution_mode == "Submit to CLI (Multi-threaded)":
                n_threads = st.slider(
                    "Number of threads",
                    min_value=1,
                    max_value=multiprocessing.cpu_count(),
                    value=default_n_threads,
                    help="Higher values = faster simulations"
                )
                st.info(f"💻 CLI will use **{n_threads}** parallel workers")
            else:
                n_threads = 1
                st.info("🌐 Browser mode uses single-threaded execution")

        # Use loaded config for simulation name and user name
        default_sim_name = ""
        default_user_name = ""
        if "loaded_config" in st.session_state and st.session_state.loaded_config:
            default_sim_name = st.session_state.loaded_config.get("simulation_name", "")
            default_user_name = st.session_state.loaded_config.get("user_name", "")

        st.markdown("")

        col_name1, col_name2 = st.columns(2)
        with col_name1:
            simulation_name = st.text_input(
                "Simulation name",
                value=default_sim_name,
                placeholder="e.g., Case_2024_001"
            )
        with col_name2:
            user_name = st.text_input(
                "Your name",
                value=default_user_name,
                help="Included in the report",
                placeholder="e.g., John Doe"
            )

        # Show a note if config was loaded
        if "loaded_config" in st.session_state and st.session_state.loaded_config:
            col_note1, col_note2 = st.columns([3, 1])
            with col_note1:
                st.info("ℹ️ Settings loaded from config.ini. You can modify them above before starting.")
            with col_note2:
                if st.button("🔄 Clear Config", help="Clear loaded configuration and use default values"):
                    st.session_state.loaded_config = None
                    st.rerun()

        st.markdown("")

        if not st.button("🚀 Start Simulation", type="primary", width='stretch'):
            return None

    # Clear loaded config after starting simulation
    if "loaded_config" in st.session_state:
        st.session_state.loaded_config = None

    input_placeholder.empty()

    if st.button("⏹️ Stop Simulation", type="secondary"):
        st.rerun()

    st.warning("⚠️ **Simulation in progress.** Do not close the browser or navigate away from this page.")

    # Create results folder
    folder_name = f"{simulation_name.replace(' ', '_').lower()}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
    results_path = Path(__file__).resolve().parent / "results" / folder_name
    results_path.mkdir(parents=True, exist_ok=True)

    # Prepare simulation parameters
    # bias_value is already set from the GUI input (None, float, or overridden by adaptive_bias)
    simulation_parameters = SimulationParameters(
        two_step_mutation_factor=two_step_mutation_factor,
        stability_window=stability_window,
        model_validity_threshold=model_validity,
        simulation_name=simulation_name,
        number_of_threads=n_threads,
        results_path=results_path,
        user_name=user_name,
        bias=bias_value,
    )

    if execution_mode == "Submit to CLI (Multi-threaded)":
        # SUBMIT TO CLI MODE
        from pedigree_lr.config_generator import generate_config_from_streamlit, run_cli_subprocess
        import json

        # Check if we have pedigree and marker set (either from upload or editors)
        if st.session_state.pedigree is None or st.session_state.marker_set is None:
            st.error("❌ Please upload or create a pedigree and haplotypes before running the simulation.")
            return None

        # Save files to temp directory
        temp_dir = results_path / "temp_files"
        temp_dir.mkdir(exist_ok=True)

        # Generate or use existing pedigree file content
        if st.session_state.pedigree_file_content is None:
            # Pedigree was created via Pedigree Builder - generate TGF content
            pedigree_content = st.session_state.pedigree.write_to_tgf().decode("utf-8")
            pedigree_extension = ".tgf"
        else:
            # Pedigree was uploaded - use existing content
            pedigree_content = st.session_state.pedigree_file_content
            pedigree_extension = st.session_state.pedigree_file_extension

        # Save pedigree file
        pedigree_path = temp_dir / f"pedigree{pedigree_extension}"
        with open(pedigree_path, "w") as f:
            f.write(pedigree_content)

        # Generate or use existing haplotypes file content
        if st.session_state.haplotypes_file_content is None:
            # Haplotypes were created via Haplotype Editor - generate JSON content
            # Build haplotypes dict from pedigree individuals
            haplotypes_dict = {}

            # Add TRACE if it exists
            if st.session_state.get("trace", None) is not None:
                trace_haplotype = st.session_state.trace
                haplotypes_dict["TRACE"] = {}
                for marker_name, alleles in trace_haplotype.alleles.items():
                    allele_strs = [str(a) for a in alleles]
                    haplotypes_dict["TRACE"][marker_name] = ";".join(allele_strs)

            # Add known individuals
            for individual in st.session_state.pedigree.get_known_individuals():
                haplotypes_dict[individual.name] = {}
                for marker_name, alleles in individual.haplotype.alleles.items():
                    allele_strs = [str(a) for a in alleles]
                    haplotypes_dict[individual.name][marker_name] = ";".join(allele_strs)

            # Add suspect if in standard mode
            if not trace_mode and st.session_state.get("suspect", None) is not None:
                suspect = st.session_state.pedigree.get_suspect()
                if suspect and suspect.name not in haplotypes_dict:
                    haplotypes_dict[suspect.name] = {}
                    for marker_name, alleles in suspect.haplotype.alleles.items():
                        allele_strs = [str(a) for a in alleles]
                        haplotypes_dict[suspect.name][marker_name] = ";".join(allele_strs)

            haplotypes_content = json.dumps(haplotypes_dict, indent=4)
        else:
            # Haplotypes were uploaded - use existing content
            haplotypes_content = st.session_state.haplotypes_file_content

        # Save haplotypes file
        haplotypes_path = temp_dir / "haplotypes.json"
        with open(haplotypes_path, "w") as f:
            f.write(haplotypes_content)

        # Store results_path in session state so we know where to find the PDF
        st.session_state.cli_results_path = results_path

        # Generate config.ini
        config_path = generate_config_from_streamlit(
            simulation_name=simulation_name,
            user_name=user_name,
            pedigree_file_path=str(pedigree_path),
            haplotypes_file_path=str(haplotypes_path),
            marker_set=st.session_state.marker_set,
            suspect_name=None if trace_mode else st.session_state.suspect,
            excluded_individuals=excluded_individuals,
            simulation_parameters=simulation_parameters,
            output_dir=results_path,
            trace_mode=trace_mode,
        )

        # Run CLI in subprocess
        progress_container = st.empty()
        progress_container.info("🚀 Launching CLI subprocess with multiprocessing...")

        # Accumulate CLI output in a list (thread-safe since only appending)
        cli_output_lines = []

        def update_progress(line):
            """Callback to accumulate CLI output (called from background thread)"""
            # Just accumulate - don't call Streamlit widgets from threads!
            cli_output_lines.append(line)

        return_code, stdout, stderr = run_cli_subprocess(
            config_path=config_path,
            skip_inside=skip_inside,
            skip_outside=skip_outside,
            trace_mode=trace_mode,
            adaptive_bias=adaptive_bias,
            progress_callback=update_progress,
        )

        # Now display the output after subprocess completes
        st.text_area("CLI Output:", value="\n".join(cli_output_lines), height=400, disabled=True)

        if return_code == 0:
            progress_container.success("✅ **Simulation completed successfully!**")

            # Get the actual results path (stored in session state before running CLI)
            actual_results_path = st.session_state.get("cli_results_path", results_path)

            # Load results from files - look for report PDFs
            # The CLI saves: simulationname_report.pdf or simulationname_trace_report.pdf
            if trace_mode:
                report_filename = f"{simulation_name}_trace_report.pdf"
            else:
                report_filename = f"{simulation_name}_report.pdf"

            pdf_path = actual_results_path / report_filename

            # Try multiple approaches to find the PDF (WSL compatibility)
            pdf_data = None
            found_path = None

            # Method 1: Direct path check
            if pdf_path.exists():
                try:
                    with open(pdf_path, "rb") as f:
                        pdf_data = f.read()
                    found_path = pdf_path
                except Exception as e:
                    st.warning(f"Could not read PDF at expected location: {e}")

            # Method 2: Search for any report PDF in results directory
            if pdf_data is None:
                try:
                    pdf_files = list(actual_results_path.glob("*report.pdf"))
                    if pdf_files:
                        # Prefer the one matching our simulation name
                        for pdf_file in pdf_files:
                            if simulation_name in pdf_file.name:
                                found_path = pdf_file
                                break
                        # Otherwise take the first one
                        if found_path is None:
                            found_path = pdf_files[0]

                        with open(found_path, "rb") as f:
                            pdf_data = f.read()
                except Exception as e:
                    st.warning(f"Could not locate PDF in results directory: {e}")

            # Method 3: List directory contents for debugging
            if pdf_data is None:
                try:
                    all_files = list(actual_results_path.glob("*.pdf"))
                    if all_files:
                        st.info(f"Found {len(all_files)} PDF file(s) in results directory:")
                        for f in all_files:
                            st.text(f"  - {f.name}")
                        # Try to read the first PDF
                        with open(all_files[0], "rb") as f:
                            pdf_data = f.read()
                        found_path = all_files[0]
                    else:
                        st.error(f"No PDF files found in: {actual_results_path}")
                        st.text(f"Looking in: {actual_results_path}")
                        st.text(f"Directory exists: {actual_results_path.exists()}")
                        if actual_results_path.exists():
                            st.text(f"Directory contents: {list(actual_results_path.glob('*'))}")
                except Exception as e:
                    st.error(f"Error searching for PDFs: {e}")

            # Display download button if we found the PDF
            if pdf_data is not None and found_path is not None:
                st.download_button(
                    label="Download simulation report",
                    data=pdf_data,
                    file_name=found_path.name,
                    mime="application/pdf",
                    help="Download the report as a PDF file."
                )
            else:
                st.error("❌ No PDF report was generated or could not be located.")
        else:
            progress_container.error(f"❌ **Simulation failed** with return code {return_code}")
            st.error("**Error output:**")
            st.code(stderr)

        # Show full output
        with st.expander("📄 Full CLI Output"):
            st.code(stdout)

        if st.button("🔄 New Simulation", type="primary"):
            st.rerun()

        return None

    else:
        # RUN IN BROWSER MODE (existing code)
        progress_placeholder = st.empty()
        log_placeholder = st.empty()

        reporter = StreamlitReporter(
            log_container=log_placeholder.container(height=600, border=True),
            progress_container=progress_placeholder.container()
        )

        simulation_result = run_simulation(
            input_pedigree=st.session_state.pedigree,
            suspect_name=None if trace_mode else st.session_state.suspect,
            marker_set=st.session_state.marker_set,
            simulation_parameters=simulation_parameters,
            reporter=reporter,
            skip_inside=skip_inside,
            skip_outside=skip_outside,
            trace_mode=trace_mode,
            trace=st.session_state.get("trace", None) if trace_mode else None,
            adaptive_bias=adaptive_bias,
        )

        progress_placeholder.empty()

        st.text("Results:")

        # Display results based on mode
        if trace_mode:
            # Trace mode: Show normalized per-individual probabilities
            st.subheader("Trace Donor Identification Results")

            # Normalize and display probabilities (including outside match probability)
            from pedigree_lr.reporting import normalize_probabilities
            normalized_probs = normalize_probabilities(
                simulation_result.per_individual_probabilities,
                simulation_result.outside_match_probability
            )

            if normalized_probs:
                st.success(f"Most Likely Donor: **{normalized_probs[0][0]}** with {normalized_probs[0][1]:.4f}% probability")

                # Create dataframe for display
                trace_df = pd.DataFrame(normalized_probs, columns=['Individual', 'Normalized Match Probability (%)'])
                trace_df['Rank'] = range(1, len(trace_df) + 1)
                trace_df = trace_df[['Rank', 'Individual', 'Normalized Match Probability (%)']]
                st.dataframe(trace_df, width='stretch')
            else:
                st.warning("No per-individual probabilities available for trace mode.")
        else:
            # Standard mode: Show full results
            try:
                proposal_distribution_dataframe = pd.DataFrame(list(simulation_result.inside_match_probability.items()),
                                                               columns=['Number of matches', 'Probability'])
                proposal_distribution_dataframe.set_index('Number of matches', inplace=True)
                st.dataframe(proposal_distribution_dataframe.sort_values(by='Number of matches'))
            except Exception as e:
                st.error(f"Error displaying proposal distribution: {e}")

            st.text("Outside match probability:")
            st.write(f"{simulation_result.outside_match_probability:.4f}")

        # Generate appropriate report based on mode
        if trace_mode:
            report_bytes = create_trace_mode_report(
                simulation_result,
                trace=st.session_state.get("trace", None)
            )
            report_label = "Download trace donor identification report"
            report_filename = f"{simulation_result.simulation_parameters.simulation_name}_trace_report.pdf"
        else:
            report_bytes = create_html_pdf_report(simulation_result)
            report_label = "Download simulation report"
            report_filename = f"{simulation_result.simulation_parameters.simulation_name}_report.pdf"

        # write report to results directory
        report_path = simulation_result.simulation_parameters.results_path / report_filename
        with open(report_path, "wb") as f:
            f.write(report_bytes)

        st.download_button(
            label=f"📥 {report_label}",
            data=report_bytes,
            file_name=report_filename,
            mime="application/pdf",
            help="Download the report as a PDF file.",
            on_click="ignore"
        )

        if st.button("🔄 New Simulation", type="primary"):
            st.rerun()

        return simulation_result


if __name__ == '__main__':
    if "marker_set" not in st.session_state:
        st.session_state.marker_set = None
    if "pedigree" not in st.session_state:
        st.session_state.pedigree = None
    if "trace" not in st.session_state:
        st.session_state.trace = None
    if "pedigree_file_content" not in st.session_state:
        st.session_state.pedigree_file_content = None
    if "pedigree_file_extension" not in st.session_state:
        st.session_state.pedigree_file_extension = None
    if "haplotypes_file_content" not in st.session_state:
        st.session_state.haplotypes_file_content = None
    if "selected_marker_set" not in st.session_state:
        st.session_state.selected_marker_set = None

    with st.sidebar:
        st.markdown("### 📁 Upload Files")

        selected_marker_set = st.selectbox(
            "Marker set",
            get_marker_set_names(),
            index=0,
            help="Select the marker set defining markers and their mutation rates.",
        )

        pedigree_file = st.file_uploader(
            "Pedigree file",
            type=["tgf", "ped"],
            help="Upload TGF or PED format. Node labels must match haplotypes file names.",
            accept_multiple_files=False
        )

        haplotypes_file = st.file_uploader(
            "Haplotypes file",
            type=["json"],
            help="Upload JSON format. Create using the Haplotype Editor.",
            accept_multiple_files=False
        )

        st.divider()

        st.markdown("### 📂 Or Load Previous Simulation")
        config_file = st.file_uploader(
            "Config.ini file",
            type=["ini"],
            help="Restore all settings and files from a previous simulation.",
            accept_multiple_files=False,
            key="config_uploader"
        )

        if config_file is not None:
            if st.button("Load simulation from config", type="primary", key="load_config_button"):
                try:
                    # Parse the config file
                    from configparser import ConfigParser
                    config = ConfigParser()
                    config.optionxform = str  # Preserve case
                    config.read_string(config_file.getvalue().decode("utf-8"))

                    # Extract paths from config
                    pedigree_path = Path(config["pedigree"]["path"])
                    haplotypes_path = Path(config["pedigree"]["known_haplotypes"])
                    marker_set_path = Path(config["pedigree"]["marker_set"])

                    # Check if files exist
                    missing_files = []
                    if not pedigree_path.exists():
                        missing_files.append(f"Pedigree: {pedigree_path}")
                    if not haplotypes_path.exists():
                        missing_files.append(f"Haplotypes: {haplotypes_path}")
                    if not marker_set_path.exists():
                        missing_files.append(f"Marker set: {marker_set_path}")

                    if missing_files:
                        st.error("❌ **Cannot load simulation** - the following files are missing:")
                        for missing in missing_files:
                            st.error(f"  • {missing}")
                        st.info("ℹ️ Ensure the config.ini file paths are correct and files are accessible.")
                    else:
                        # Load marker set
                        from pedigree_lr.models import MarkerSet
                        marker_set = MarkerSet()
                        with open(marker_set_path, "r") as f:
                            marker_set.read_marker_set_from_file(f)
                        st.session_state.marker_set = marker_set

                        # Determine marker set name (try to match with database)
                        marker_set_names = get_marker_set_names()
                        matched_marker_set_name = None
                        for ms_name in marker_set_names:
                            test_ms = load_marker_set_from_database(ms_name)
                            if len(test_ms.markers) == len(marker_set.markers):
                                if all(m1.name == m2.name for m1, m2 in zip(test_ms.markers, marker_set.markers)):
                                    matched_marker_set_name = ms_name
                                    break

                        if matched_marker_set_name:
                            st.session_state.selected_marker_set = matched_marker_set_name
                        else:
                            st.session_state.selected_marker_set = marker_set_names[0] if marker_set_names else None

                        # Load pedigree
                        pedigree_extension = pedigree_path.suffix
                        with open(pedigree_path, "r") as f:
                            pedigree_content = f.read()
                            stringio = StringIO(pedigree_content)
                            st.session_state.pedigree = load_pedigree_from_upload(stringio, pedigree_extension)
                            st.session_state.pedigree_file_content = pedigree_content
                            st.session_state.pedigree_file_extension = pedigree_extension

                        # Load haplotypes
                        with open(haplotypes_path, "r") as f:
                            haplotypes_content = f.read()
                            stringio = StringIO(haplotypes_content)
                            trace_haplotype = st.session_state.pedigree.read_known_haplotypes_from_file(
                                stringio, st.session_state.marker_set
                            )
                            st.session_state.haplotypes_file_content = haplotypes_content

                            # Store trace if present
                            if trace_haplotype is not None:
                                st.session_state.trace = trace_haplotype
                            else:
                                st.session_state.trace = None

                        # Store config parameters in session state for later use
                        st.session_state.loaded_config = {
                            "simulation_name": config["pedigree"].get("simulation_name", ""),
                            "user_name": config["pedigree"].get("user_name", ""),
                            "suspect": config["pedigree"].get("suspect", None),
                            "exclude_individuals": config["pedigree"].get("exclude_individuals", "").split(",") if config["pedigree"].get("exclude_individuals", "") else [],
                            "two_step_mutation_fraction": float(config["pedigree"]["two_step_mutation_fraction"]),
                            "batch_length": int(config["pedigree"]["batch_length"]),
                            "convergence_criterion": float(config["pedigree"]["convergence_criterion"]),
                            "number_of_threads": int(config["pedigree"]["number_of_threads"]),
                            "bias": float(config["pedigree"]["bias"]) if "bias" in config["pedigree"] else None,
                            "trace_mode": trace_haplotype is not None,
                        }

                        st.success("✅ **Configuration loaded successfully!**")
                        st.info("ℹ️ All settings and files have been restored. Review the settings below and start the simulation.")
                        st.rerun()

                except Exception as e:
                    st.error(f"❌ **Error loading configuration:** {str(e)}")
                    import traceback
                    st.code(traceback.format_exc())

        st.divider()

        if st.button("📤 Upload Files", type="primary", width='stretch'):
            if selected_marker_set:
                st.session_state.marker_set = load_marker_set_from_database(selected_marker_set)
                st.session_state.selected_marker_set = selected_marker_set

            if pedigree_file is not None:
                file_extension = Path(pedigree_file.name).suffix
                stringio = StringIO(pedigree_file.getvalue().decode("utf-8"))
                st.session_state.pedigree = load_pedigree_from_upload(stringio, file_extension)

                # Store file content for CLI mode
                st.session_state.pedigree_file_content = pedigree_file.getvalue().decode("utf-8")
                st.session_state.pedigree_file_extension = file_extension

                if st.session_state.pedigree is not None:
                    st.success("✅ Pedigree uploaded")

            if haplotypes_file is not None:
                stringio = StringIO(haplotypes_file.getvalue().decode("utf-8"))
                trace_haplotype = st.session_state.pedigree.read_known_haplotypes_from_file(
                    stringio, st.session_state.marker_set
                )

                # Store file content for CLI mode
                st.session_state.haplotypes_file_content = haplotypes_file.getvalue().decode("utf-8")

                # Store trace in session state
                if trace_haplotype is not None:
                    st.session_state.trace = trace_haplotype
                    st.success("✅ Haplotypes uploaded (with TRACE)")
                else:
                    st.session_state.trace = None
                    st.success("✅ Haplotypes uploaded")

        st.divider()

        st.markdown("### 📥 Example Files")
        with open(r"examples/example.tgf") as file:
            st.download_button(
                "Download example pedigree",
                file,
                "example_pedigree.tgf",
                help="Download an example pedigree file in TGF format",
                width='stretch'
            )

    # Main content area
    st.title("🧬 MatchY - Y-STR Pedigree Likelihood Ratio Analysis")
    st.markdown("Perform likelihood ratio analysis for Y-STR pedigree-based forensic casework.")
    st.markdown("---")

    if st.session_state.pedigree is not None:
        result = render_simulation()
    else:
        with st.container():
            st.info("### 👋 Welcome to MatchY!")
            st.markdown("""
            To get started, you need to provide:
            1. **Pedigree file** - The family structure
            2. **Haplotypes file** - Y-STR profiles for known individuals

            You can:
            - 📤 **Upload files** using the sidebar
            - 🏗️ **Build a pedigree** using the Pedigree Builder page
            - 📝 **Create haplotypes** using the Haplotypes Editor page
            """)

            col_help1, col_help2, col_help3 = st.columns(3)
            with col_help1:
                st.markdown("##### 📚 Need Help?")
                st.markdown("Check the example files in the sidebar")
            with col_help2:
                st.markdown("##### 🎯 Quick Start")
                st.markdown("Upload a pedigree and haplotypes file")
            with col_help3:
                st.markdown("##### 🔬 Trace Mode")
                st.markdown("Include TRACE profile for donor identification")
