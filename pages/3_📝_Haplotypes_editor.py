from io import StringIO
import streamlit as st
from pedigree_lr.data import get_marker_set_names, load_marker_set_from_database
import pandas as pd
from pathlib import Path
import json

st.set_page_config(
    page_title="Haplotypes editor",
    page_icon="📝",
    layout="wide",
    initial_sidebar_state="expanded",
)

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

st.title("🧬 Haplotype Editor")
st.markdown("Create and edit Y-STR haplotype profiles for individuals in your pedigree.")
st.markdown("---")

# Configuration Section
with st.container():
    st.subheader("⚙️ Configuration")
    col_config1, col_config2 = st.columns([2, 3])
    with col_config1:
        selected_marker_set = st.selectbox("Select marker set",
                                           get_marker_set_names(),
                                           index=0,
                                           help="Select the marker set you want to use for editing haplotypes.",
                                           )
    with col_config2:
        with st.expander("ℹ️ Allele Format Guidelines"):
            st.markdown("""
            **Allowed allele formats:**
            - Single allele: `x` (e.g., `15`)
            - Decimal allele: `x.y` (e.g., `15.2`)
            - Multi-copy: `x;y` (e.g., `15;16`)
            - Multi-copy with decimal: `x.y;z.w` (e.g., `15.2;16.3`)

            Where x, y, z, w are integers and `;` separates multi-copy alleles.
            """)

haplotype_markers = load_marker_set_from_database(selected_marker_set)
haplotype_markers_df = pd.DataFrame([marker.name for marker in haplotype_markers.markers],
                                    columns=["Markers"])

# Initialize or update haplotype_markers_df in session state
if "haplotype_markers_df" not in st.session_state:
    st.session_state.haplotype_markers_df = haplotype_markers_df
else:
    # If marker set changed, reset the dataframe
    if st.session_state.haplotype_markers_df.shape[0] != haplotype_markers_df.shape[0]:
        st.session_state.haplotype_markers_df = haplotype_markers_df

# Feature 1: Load known haplotypes from pedigree object if available
# Check if pedigree exists and if the selected marker set matches the one used by the pedigree
if "pedigree" in st.session_state and st.session_state.pedigree is not None:
    # Check if we have a marker set stored and if it matches the selected one
    if "selected_marker_set" in st.session_state and st.session_state.selected_marker_set == selected_marker_set:
        # Check if we need to populate the dataframe with known haplotypes from pedigree
        # Only populate if the dataframe only has the Markers column (no individuals added yet)
        if st.session_state.haplotype_markers_df.shape[1] == 1:
            # Get known individuals from pedigree
            known_individuals = st.session_state.pedigree.get_known_individuals()
            suspect = st.session_state.pedigree.get_suspect()

            # Combine known individuals and suspect
            all_known = known_individuals.copy()
            if suspect is not None:
                all_known.append(suspect)

            # Populate dataframe with known haplotypes
            if len(all_known) > 0:
                for individual in all_known:
                    # Extract haplotype for this individual
                    haplotype_values = []
                    for marker_name in st.session_state.haplotype_markers_df["Markers"]:
                        if marker_name in individual.haplotype.alleles:
                            alleles = individual.haplotype.alleles[marker_name]
                            # Convert alleles to string format (e.g., "15;16" for multi-copy)
                            allele_strs = [str(a) for a in alleles]
                            haplotype_values.append(";".join(allele_strs))
                        else:
                            haplotype_values.append("")

                    # Add column for this individual
                    st.session_state.haplotype_markers_df[individual.name] = haplotype_values

            # Check if TRACE exists in session state and add it
            if "trace" in st.session_state and st.session_state.trace is not None:
                trace_haplotype = st.session_state.trace
                trace_values = []
                for marker_name in st.session_state.haplotype_markers_df["Markers"]:
                    if marker_name in trace_haplotype.alleles:
                        alleles = trace_haplotype.alleles[marker_name]
                        allele_strs = [str(a) for a in alleles]
                        trace_values.append(";".join(allele_strs))
                    else:
                        trace_values.append("")

                # Add TRACE as the first column after Markers (if not already present)
                if "TRACE" not in st.session_state.haplotype_markers_df.columns:
                    st.session_state.haplotype_markers_df.insert(1, "TRACE", trace_values)

st.markdown("")  # Spacing

# Haplotype Data Section
with st.container():
    st.subheader("📊 Haplotype Data")

    # Check if TRACE exists and apply custom styling
    if "TRACE" in st.session_state.haplotype_markers_df.columns:
        st.markdown("""
        <style>
        /* Style for TRACE column header */
        [data-testid="stDataFrameResizable"] th:nth-child(2) {
            background-color: #dc2626 !important;
            color: white !important;
            font-weight: 700 !important;
        }
        </style>
        """, unsafe_allow_html=True)

    # Render data editor FIRST to capture user edits
    edited_df = st.data_editor(st.session_state.haplotype_markers_df,
                               hide_index=True,
                               disabled=["Markers"],
                               key="haplotype_editor",
                               use_container_width=True
                               )

    # Validation messages
    if edited_df.isin([""]).any().any():
        st.error("⚠️ Haplotypes cannot contain empty values. Please fill in all fields or change the marker set.")
    else:
        validation_errors = []
        for index, row in edited_df.iterrows():
            copy_numbers = []
            for col in edited_df.columns[1:]:
                alleles = row[col].split(";")
                copy_numbers.append(len(alleles))
                for allele in alleles:
                    values = allele.split(".")
                    for value in values:
                        if not value.isdigit():
                            validation_errors.append(f"Invalid allele value at marker '{row['Markers']}'. Allele values must be in the form 'x.y' or 'x'")
                            break
            if len(set(copy_numbers)) > 1:
                validation_errors.append(f"Different number of copies at marker '{row['Markers']}'. All individuals must have the same number of copies.")

        if validation_errors:
            for error in validation_errors[:3]:  # Show first 3 errors
                st.error(f"⚠️ {error}")
            if len(validation_errors) > 3:
                st.error(f"... and {len(validation_errors) - 3} more validation errors.")

st.markdown("---")

# Add/Remove Profiles Section
with st.container():
    st.subheader("➕ Add or Remove Profiles")

    # Initialize file uploader key in session state to allow clearing after use
    if "haplotype_uploader_key" not in st.session_state:
        st.session_state.haplotype_uploader_key = 0

    # Create tabs for different actions
    tab1, tab2, tab3 = st.tabs(["👤 Add Individual", "🔬 Add TRACE Profile", "🗑️ Remove Profiles"])

    with tab1:
        st.markdown("##### Add a new individual to the haplotype table")

        col_upload1, col_select1 = st.columns([2, 3])

        with col_upload1:
            uploaded_haplotype_file_individual = st.file_uploader(
                "Upload haplotype file (optional)",
                type=["txt", "csv"],
                help="Upload a comma-separated file with marker and alleles columns. Leave empty to create an empty haplotype or copy from the last column.",
                key=f"haplotype_uploader_individual_{st.session_state.haplotype_uploader_key}"
            )

        with col_select1:
            if "pedigree" in st.session_state and st.session_state.pedigree is not None:
                available_individuals = [individual.name for individual in st.session_state.pedigree.individuals if individual.name not in edited_df.columns]
                new_individual = st.selectbox(
                    "Select individual from pedigree",
                    available_individuals,
                    index=0 if len(available_individuals) > 0 else None,
                    help="Select an individual from your loaded pedigree."
                )
            else:
                st.warning("⚠️ No pedigree loaded. Manual entry mode.")
                new_individual = st.text_input("Enter individual name", help="Manually enter a new individual name (not recommended).")

        if st.button("➕ Add Individual", type="primary", disabled=new_individual == "" or new_individual is None):
            if new_individual:
                if new_individual.upper() == "TRACE":
                    st.error("❌ Cannot add 'TRACE' as a regular individual. Use the 'Add TRACE Profile' tab.")
                elif new_individual in edited_df.columns:
                    st.warning(f"⚠️ Individual '{new_individual}' already exists.")
                else:
                    file_was_uploaded = uploaded_haplotype_file_individual is not None
                    if file_was_uploaded:
                        uploaded_haplotype_df = pd.read_csv(uploaded_haplotype_file_individual, sep=",", header=0, names=["marker", "alleles"])
                        haplotype_dict = dict(zip(uploaded_haplotype_df.marker, uploaded_haplotype_df.alleles))
                        new_haplotype = [haplotype_dict.get(marker, "") for marker in edited_df["Markers"]]
                        edited_df[new_individual] = new_haplotype
                    elif len(edited_df.columns) > 1:
                        edited_df[new_individual] = edited_df.iloc[:, -1].copy()
                    else:
                        edited_df[new_individual] = [""] * len(edited_df)

                    st.session_state.haplotype_markers_df = edited_df
                    if file_was_uploaded:
                        st.session_state.haplotype_uploader_key += 1
                    st.rerun()

    with tab2:
        st.markdown("##### Add a TRACE profile for trace donor identification analysis")
        st.info("💡 The TRACE profile represents an unknown DNA sample to be compared against pedigree members. If no file is uploaded, it will copy the last individual's haplotype as a starting point.")

        uploaded_haplotype_file_trace = st.file_uploader(
            "Upload TRACE haplotype file (optional)",
            type=["txt", "csv"],
            help="Upload a comma-separated file with marker and alleles columns. If not provided, the last individual's haplotype will be copied.",
            key=f"haplotype_uploader_trace_{st.session_state.haplotype_uploader_key}"
        )

        if st.button("🔬 Add TRACE Profile", type="primary", disabled="TRACE" in edited_df.columns):
            if "TRACE" not in edited_df.columns:
                file_was_uploaded = uploaded_haplotype_file_trace is not None
                if file_was_uploaded:
                    uploaded_haplotype_df = pd.read_csv(uploaded_haplotype_file_trace, sep=",", header=0, names=["marker", "alleles"])
                    haplotype_dict = dict(zip(uploaded_haplotype_df.marker, uploaded_haplotype_df.alleles))
                    trace_haplotype = [haplotype_dict.get(marker, "") for marker in edited_df["Markers"]]
                else:
                    # If there are existing individuals, copy from the last one; otherwise create empty
                    if len(edited_df.columns) > 1:
                        trace_haplotype = edited_df.iloc[:, -1].copy().tolist()
                    else:
                        trace_haplotype = [""] * len(edited_df)

                edited_df.insert(1, "TRACE", trace_haplotype)
                st.session_state.haplotype_markers_df = edited_df
                if file_was_uploaded:
                    st.session_state.haplotype_uploader_key += 1
                st.rerun()

        if "TRACE" in edited_df.columns:
            st.success("✅ TRACE profile is present in the table.")

    with tab3:
        st.markdown("##### Remove individuals or TRACE profile from the table")

        col_remove1, col_remove2 = st.columns(2)

        with col_remove1:
            if "pedigree" in st.session_state and st.session_state.pedigree is not None:
                available_to_remove = [col for col in edited_df.columns[1:] if col != "TRACE"]
                if len(available_to_remove) > 0:
                    individual_to_remove = st.selectbox("Select individual to remove", available_to_remove)
                    if st.button("🗑️ Remove Selected Individual", type="secondary"):
                        edited_df = edited_df.drop(columns=[individual_to_remove])
                        st.session_state.haplotype_markers_df = edited_df
                        st.rerun()
                else:
                    st.info("No individuals to remove.")
            else:
                st.info("Load a pedigree to remove individuals.")

        with col_remove2:
            if "TRACE" in edited_df.columns:
                if st.button("🗑️ Remove TRACE Profile", type="secondary"):
                    edited_df = edited_df.drop(columns=["TRACE"])
                    st.session_state.haplotype_markers_df = edited_df
                    st.rerun()
            else:
                st.info("No TRACE profile to remove.")

st.markdown("---")

# Actions Section
with st.container():
    st.subheader("💾 Save and Export")

    col_action1, col_action2, col_action3 = st.columns(3)

    with col_action1:
        st.markdown("**Save Changes**")
        if st.button("💾 Save Changes to Session", type="secondary", use_container_width=True):
            st.session_state.haplotype_markers_df = edited_df
            st.success("✅ Changes saved to session.")


    with col_action2:
        st.markdown("**Download Haplotypes**")

        def write_to_json(df):
            haplotypes_dict = {}
            for individual in df.columns[1:]:
                haplotypes_dict[individual] = {}

            for i, marker in df.iterrows():
                for individual in df.columns[1:]:
                    haplotypes_dict[individual][marker["Markers"]] = marker[individual]

            return json.dumps(haplotypes_dict, indent=4)

        st.download_button(
            "📥 Download as JSON",
            data=write_to_json(edited_df),
            file_name=f"haplotypes.json",
            help="Download the haplotypes as a JSON file for later use.",
            use_container_width=True
        )

    with col_action3:
        st.markdown("**Load to Simulation**")
        if st.button(
            "🚀 Load to Simulation",
            type="primary",
            help="Load haplotypes directly to the simulation. Requires a pedigree to be loaded.",
            disabled="pedigree" not in st.session_state,
            use_container_width=True
        ):
            st.session_state.haplotype_markers_df = edited_df

            # Generate JSON content from edited dataframe
            json_content = write_to_json(edited_df)
            stringio = StringIO(json_content)
            st.session_state.marker_set = load_marker_set_from_database(selected_marker_set)

            # Store the JSON content in session state for CLI mode
            st.session_state.haplotypes_file_content = json_content

            # Read haplotypes and capture TRACE if present
            trace_haplotype = st.session_state.pedigree.read_known_haplotypes_from_file(stringio, st.session_state.marker_set)

            # Store TRACE in session state if it was found
            if trace_haplotype is not None:
                st.session_state.trace = trace_haplotype
                st.success("✅ Haplotypes loaded to simulation (including TRACE profile).")
                st.info("➡️ Go to Home page to start the simulation.")
            else:
                st.session_state.trace = None
                st.success("✅ Haplotypes loaded to simulation.")
                st.info("➡️ Go to Home page to start the simulation.")