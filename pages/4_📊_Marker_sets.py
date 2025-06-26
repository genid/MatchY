import math

import streamlit as st
import pandas as pd
from pathlib import Path
import json

st.set_page_config(
    page_title="Marker sets",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded",
)

mutation_rates_path = Path(__file__).resolve().parent.parent / "data" / "mutation_rates.csv"

mutation_rates_df = pd.read_csv(mutation_rates_path, index_col=0, header=0, names=["Marker", "Mutation rate"])

col1, col2 = st.columns(2, gap="large")

col1.header("Mutation rates")

st.info("Double click on a cell to edit it. Empty rows will be removed from the table.")
edited_df = col1.data_editor(
    mutation_rates_df
)

if col1.button("Save changes",
             type="primary", ):
    edited_df = edited_df.dropna()
    edited_df.to_csv(mutation_rates_path, index=True, header=True)
    st.rerun()

col3, col4 = col1.columns(2)

new_marker_name = col3.text_input("Enter new marker name")
new_marker_rate = col4.number_input("Enter new marker mutation rate",
                                    min_value=0.0000,
                                    max_value=1.0000,
                                    value=0.0001,
                                    step=0.0001,
                                    format="%.8f")

if col1.button("Add new marker",
             type="primary", ):
    if new_marker_name and new_marker_rate:
        if new_marker_name in edited_df.index:
            st.error("Marker already exists. Please choose a different name.")
        elif new_marker_rate < 0.0000 or new_marker_rate > 1.0000:
            st.error("Mutation rate must be between 0 and 1.")
        else:
            new_row = pd.DataFrame([[new_marker_name, new_marker_rate]], columns=["Marker", "Mutation rate"])
            new_row.set_index("Marker", inplace=True)
            edited_df = pd.concat([edited_df, new_row], ignore_index=False)
            edited_df.to_csv(mutation_rates_path, index=True, header=True)
            st.rerun()
    else:
        st.error("Please enter both marker name and mutation rate.")

if col1.button("Reset mutation rates to default (cannot be undone!)", type="tertiary"):
    backup_path = Path(__file__).resolve().parent.parent / "data" / "mutation_rates_backup.csv"
    default_df = pd.read_csv(backup_path, index_col=0, header=0, names=["Marker", "Mutation rate"])
    edited_df = default_df
    edited_df.to_csv(mutation_rates_path, index=True, header=True)
    st.rerun()

col2.header("Marker sets")

kits_path = Path(__file__).resolve().parent.parent / "data" / "kits.json"

with open(kits_path, "r") as file:
    kits = json.load(file)

kit_names = list(kits.keys())
kit_names.append("Add new kit")

selected_kit = col2.selectbox("Select a kit", kit_names)

if selected_kit == "Add new kit":
    st.info("Enter a name for the new kit below and select the markers to include in it.")
    selected_kit = col2.text_input("Enter new kit name")
    kits[selected_kit] = []

selected_kit_markers = kits[selected_kit]

marker_multiselect = col2.multiselect(
    "Select markers",
    options=mutation_rates_df.index.tolist(),
    default=selected_kit_markers,
    help="Select the markers to include in the kit. "
         "The markers must be in the mutation rates table.",
)

overall_mutation_rate = 1 - math.prod([1 - mut_rate for mut_rate in mutation_rates_df.loc[marker_multiselect, "Mutation rate"]])
col2.info(f"Overall mutation rate for selected markers: {overall_mutation_rate:.4f} mutation(s) per generation (in at least one marker).")

if col2.button("Save kit",
             type="primary", ):
    if selected_kit and marker_multiselect and len(marker_multiselect) > 0:
        kits[selected_kit] = marker_multiselect
        with open(kits_path, "w") as file:
            json.dump(kits, file, indent=4)
        st.rerun()
    else:
        st.error("Please select a kit name and at least one marker.")

if col2.button("Delete kit"):
    if selected_kit and selected_kit in kits:
        del kits[selected_kit]
        with open(kits_path, "w") as file:
            json.dump(kits, file, indent=4)
        st.rerun()
    else:
        st.error("Please select a kit to delete.")

if col2.button("Reset marker sets to default (cannot be undone!)", type="tertiary"):
    backup_path = Path(__file__).resolve().parent.parent / "data" / "kits_backup.json"
    with open(backup_path, "r") as file:
        default_kits = json.load(file)
    kits = default_kits
    with open(kits_path, "w") as file:
        json.dump(kits, file, indent=4)
    st.rerun()