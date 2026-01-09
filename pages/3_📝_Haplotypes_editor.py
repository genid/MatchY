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

selected_marker_set = st.selectbox("Select marker set",
                                   get_marker_set_names(),
                                   index=0,
                                   help="Select the marker set you want to use for editing haplotypes.",
                                   )

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

st.warning(
    "Allowed allele values are: 'x', 'x.y', 'x.y;z', 'x.y;z.w', where x, y, z, w are integers and ; is the seperator for multi-copy alleles.")

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
                           key="haplotype_editor"
                           )

# Add TRACE button after the editor so we can preserve edits
st.info("You can add individuals or a TRACE profile to the haplotypes table. Upload a CSV file (optional) to populate the haplotype values automatically.")

# Initialize file uploader key in session state to allow clearing after use
if "haplotype_uploader_key" not in st.session_state:
    st.session_state.haplotype_uploader_key = 0

uploaded_haplotype_file = st.file_uploader("Upload haplotype file (or leave empty for an empty haplotype)",
                                           type=["txt", "csv"],
                                           help="Upload a comma-separated file with the haplotypes. The file should contain a header with two columns: 'marker' and 'alleles'. When no file is uploaded, an empty haplotype will be created.",
                                           key=f"haplotype_uploader_{st.session_state.haplotype_uploader_key}"
                                           )

col_trace1, col_trace2 = st.columns([1, 4])
with col_trace1:
    if st.button("Add TRACE Profile",
                 type="secondary",
                 help="Add a TRACE profile for trace mode analysis. The TRACE profile represents an unknown DNA sample to be compared against pedigree members.",
                 disabled="TRACE" in edited_df.columns):
        if "TRACE" not in edited_df.columns:
            # Parse uploaded file if present
            file_was_uploaded = uploaded_haplotype_file is not None
            if file_was_uploaded:
                uploaded_haplotype_df = pd.read_csv(uploaded_haplotype_file, sep=",", header=0,
                                                    names=["marker", "alleles"])
                haplotype_dict = dict(zip(uploaded_haplotype_df.marker, uploaded_haplotype_df.alleles))
                trace_haplotype = []
                for marker in edited_df["Markers"]:
                    if marker in haplotype_dict:
                        trace_haplotype.append(haplotype_dict[marker])
                    else:
                        trace_haplotype.append("")
            else:
                # Create empty TRACE column
                trace_haplotype = [""] * len(edited_df)

            # Add TRACE as the first column after Markers, preserving existing data from edited_df
            edited_df.insert(1, "TRACE", trace_haplotype)
            st.session_state.haplotype_markers_df = edited_df

            # Clear the file uploader if a file was used
            if file_was_uploaded:
                st.session_state.haplotype_uploader_key += 1

            st.rerun()
        else:
            st.warning("TRACE profile already exists.")

with col_trace2:
    if "TRACE" in edited_df.columns:
        st.info("TRACE profile added.")

if edited_df.isin([""]).any().any():
    st.error("Haplotypes cannot contain empty values. Please fill in all fields or change the marker set.")
else:
    for index, row in edited_df.iterrows():
        copy_numbers = []
        for col in edited_df.columns[1:]:
            alleles = row[col].split(";")
            copy_numbers.append(len(alleles))
            for allele in alleles:
                values = allele.split(".")
                for value in values:
                    if not value.isdigit():
                        st.error(f"Invalid allele value. Allele value must be in the form of 'x.y' or 'x'")
                        break
        if len(set(copy_numbers)) > 1:
            st.error(f"Different number of copies. All individuals must have the same number of copies.")
            break

st.info("You can add new individuals to the haplotypes table.")

if "pedigree" in st.session_state and st.session_state.pedigree is not None:
    # Filter out TRACE from the individual list
    available_individuals = [individual.name for individual in st.session_state.pedigree.individuals if individual.name not in edited_df.columns]
    new_individual = st.selectbox("Select individual name",
                                  available_individuals,
                                  index=0,
                                  help="Select the individual you want to add a haplotype for.",
                                  )

else:
    st.error(
        "You don't have a pedigree loaded yet. Manually entering individuals is not recommended. Use this only to create haplotype files for future usage.")
    new_individual = st.text_input("Manually enter a new individual name (do not use 'TRACE' - use the TRACE button above).")

col1, col2, col3, col4, col5, col6 = st.columns(6)

if col1.button("Add selected individual",
               type="primary",
               disabled=new_individual == "" or new_individual is None,
               help="Add the selected individual to the haplotypes table. If a file is uploaded above, it will be used to populate the haplotype values. Otherwise, an empty haplotype will be created."):
    if new_individual:
        # Prevent adding TRACE manually
        if new_individual.upper() == "TRACE":
            st.error("Cannot add 'TRACE' as a regular individual. Use the 'Add TRACE Profile' button above instead.")
        elif new_individual in edited_df.columns:
            st.warning(f"Individual {new_individual} already exists.")
        else:
            file_was_uploaded = uploaded_haplotype_file is not None
            if file_was_uploaded:
                uploaded_haplotype_df = pd.read_csv(uploaded_haplotype_file, sep=",", header=0,
                                                    names=["marker", "alleles"])
                haplotype_dict = dict(zip(uploaded_haplotype_df.marker, uploaded_haplotype_df.alleles))
                new_haplotype = []
                for marker in edited_df["Markers"]:
                    if marker in haplotype_dict:
                        new_haplotype.append(haplotype_dict[marker])
                    else:
                        new_haplotype.append("")
                edited_df[new_individual] = new_haplotype
            elif len(edited_df.columns) > 1:
                edited_df[new_individual] = edited_df.iloc[:, -1].copy()
            else:
                # Create empty haplotype
                edited_df[new_individual] = [""] * len(edited_df)

            st.session_state.haplotype_markers_df = edited_df

            # Clear the file uploader if a file was used
            if file_was_uploaded:
                st.session_state.haplotype_uploader_key += 1

            st.rerun()
    else:
        st.warning("Please enter a name for the new individual.")

if col2.button("Save changes",
               type="secondary", ):
    st.session_state.haplotype_markers_df = edited_df
    st.success("Changes saved.")


def write_to_json(df):
    haplotypes_dict = {}
    for individual in df.columns[1:]:
        haplotypes_dict[individual] = {}

    for i, marker in df.iterrows():
        for individual in df.columns[1:]:
            haplotypes_dict[individual][marker["Markers"]] = marker[individual]

    return json.dumps(haplotypes_dict, indent=4)


if col3.download_button("Download haplotypes",
                        data=write_to_json(edited_df),
                        file_name=f"haplotypes.json",
                        help="Download the haplotypes as a JSON file."):
    st.success(f"Haplotypes downloaded as haplotypes.json")

if col4.button("Load haplotypes directly to simulation",
               type="primary",
               help="Make sure that a pedigree file is loaded first.",
               disabled="pedigree" not in st.session_state):
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
        st.success("Haplotypes loaded to simulation (including TRACE profile). Go to Home to start the simulation.")
    else:
        st.session_state.trace = None
        st.success("Haplotypes loaded to simulation. Go to Home to start the simulation.")

    st.info("Make sure to download the haplotypes file if you want to use it for other simulations.")

if col5.button("Remove selected individual",
               type="secondary",
               disabled=new_individual == "" or new_individual is None or new_individual not in edited_df.columns,
               help="The currently selected individual will be removed from the table."):
    if new_individual in edited_df.columns:
        if new_individual == "TRACE":
            st.error("Cannot remove TRACE using this button. Use the '🗑️ Remove TRACE' button instead.")
        else:
            edited_df = edited_df.drop(columns=[new_individual])
            st.session_state.haplotype_markers_df = edited_df
            st.rerun()
    else:
        st.warning(f"Individual {new_individual} does not exist.")

if col6.button("Remove TRACE",
               type="secondary",
               disabled="TRACE" not in edited_df.columns,
               help="Remove the TRACE profile from the haplotypes table."):
    if "TRACE" in edited_df.columns:
        edited_df = edited_df.drop(columns=["TRACE"])
        st.session_state.haplotype_markers_df = edited_df
        st.rerun()
    else:
        st.warning("TRACE profile does not exist.")