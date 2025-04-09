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

st.info(
    "Allowed allele values are: 'x', 'x.y', 'x.y;z', 'x.y;z.w', where x, y, z, w are integers and ; is the seperator for multi-copy alleles.")

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

new_individual = st.text_input("Enter new individual name (haplotype of the last added individual will be copied).")

if st.button("Add individual",
             type="primary", ):
    if new_individual:
        if new_individual in edited_df.columns:
            st.warning(f"Individual {new_individual} already exists.")
        else:
            if len(edited_df.columns) > 1:
                edited_df[new_individual] = edited_df.iloc[:, -1].copy()
            else:
                edited_df[new_individual] = [""] * len(edited_df)
            st.session_state.haplotype_markers_df = edited_df
            st.rerun()
    else:
        st.warning("Please enter a name for the new individual.")


def write_to_json(edited_df):
    haplotypes_dict = {}
    for col in edited_df.columns[1:]:
        haplotypes_dict[col] = {}

    for index, row in edited_df.iterrows():
        for col in edited_df.columns[1:]:
            haplotypes_dict[col][row["Markers"]] = row[col]

    return json.dumps(haplotypes_dict, indent=4)

if st.download_button("Download haplotypes",
                      data=write_to_json(edited_df),
                      file_name=f"haplotypes.json",
                      help="Download the haplotypes as a JSON file."):
    st.success(f"Haplotypes downloaded as haplotypes.json")
