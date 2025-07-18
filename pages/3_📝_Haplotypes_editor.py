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

selected_marker_set = st.selectbox("Select marker set",
                                   get_marker_set_names(),
                                   index=0,
                                   help="Select the marker set you want to use for editing haplotypes.",
                                   )

haplotype_markers = load_marker_set_from_database(selected_marker_set)
haplotype_markers_df = pd.DataFrame([marker.name for marker in haplotype_markers.markers],
                                    columns=["Markers"])

if "haplotype_markers_df" not in st.session_state:
    st.session_state.haplotype_markers_df = haplotype_markers_df
else:
    if st.session_state.haplotype_markers_df.shape[0] != haplotype_markers_df.shape[0]:
        st.session_state.haplotype_markers_df = haplotype_markers_df

col1, col2 = st.columns(2)
col1.warning(
    "Allowed allele values are: 'x', 'x.y', 'x.y;z', 'x.y;z.w', where x, y, z, w are integers and ; is the seperator for multi-copy alleles.")
col2.info("You can add new individuals to the haplotypes table.")

edited_df = st.data_editor(st.session_state.haplotype_markers_df,
                           hide_index=True,
                           disabled=["Markers"],
                           )

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

if "pedigree" in st.session_state and st.session_state.pedigree is not None:
    st.info("Since you have a pedigree loaded, you can add haplotypes to these individuals.")
    new_individual = st.selectbox("Select individual name",
                                  [individual.name for individual in st.session_state.pedigree.individuals],
                                  index=0,
                                  help="Select the individual you want to add a haplotype for.",
                                  )

else:
    st.error(
        "You don't have a pedigree loaded yet. Manually entering individuals is not recommended. Use this only to create haplotype files for future usage.")
    new_individual = st.text_input("Manually enter a new individual name.")

uploaded_haplotype_file = st.file_uploader("Upload haplotype file (or leave empty to copy the last added individual)",
                                           type=["txt", "csv"],
                                           help="Upload a comma-separated file with the haplotypes for the new individual. The file should contain two columns: 'marker' and 'alleles'. When no file is uploaded, the haplotype of the last added individual will be copied.",
                                           )

col1, col2, col3, col4 = st.columns(4)

if col1.button("Add individual",
               type="primary",
               disabled=new_individual == "" or new_individual is None,
               help="Enter the name of a new individual to add to the haplotypes table."):
    if new_individual:
        if new_individual in edited_df.columns:
            st.warning(f"Individual {new_individual} already exists.")
        else:
            if uploaded_haplotype_file is not None:
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
                edited_df[new_individual] = [""] * len(edited_df)

            st.session_state.haplotype_markers_df = edited_df
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
    stringio = StringIO(write_to_json(edited_df))
    st.session_state.marker_set = load_marker_set_from_database(selected_marker_set)
    st.session_state.pedigree.read_known_haplotypes_from_file(stringio, st.session_state.marker_set)
    st.success("Haplotypes loaded to simulation. Go to Home to start the simulation.")
    st.info("Make sure to download the haplotypes file if you want to use it for other simulations.")
