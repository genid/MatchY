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

st.title("📊 Marker Sets & Mutation Rates")
st.markdown("Manage Y-STR markers, mutation rates, and create custom marker sets for your analyses.")
st.markdown("---")

# Load data
mutation_rates_path = Path(__file__).resolve().parent.parent / "data" / "mutation_rates.csv"
mutation_rates_df = pd.read_csv(mutation_rates_path, index_col=0, header=0, names=["Marker", "Mutation rate"])

kits_path = Path(__file__).resolve().parent.parent / "data" / "kits.json"
with open(kits_path, "r") as file:
    kits = json.load(file)

# Create two-column layout
col1, col2 = st.columns(2, gap="large")

# ==================== LEFT COLUMN: MUTATION RATES ====================
with col1:
    st.subheader("🧬 Mutation Rates")

    with st.container():
        st.markdown("##### Edit Mutation Rates")
        st.info("💡 Double-click on a cell to edit. Empty rows will be automatically removed.")

        edited_df = st.data_editor(
            mutation_rates_df,
            width='stretch',
            height=400
        )

        col_save, col_reset = st.columns([1, 1])

        with col_save:
            if st.button("💾 Save Changes", type="primary", width='stretch'):
                edited_df = edited_df.dropna()
                edited_df.to_csv(mutation_rates_path, index=True, header=True)
                st.success("✅ Mutation rates saved!")
                st.rerun()

        with col_reset:
            if st.button("🔄 Reset to Default", type="secondary", width='stretch', help="⚠️ This action cannot be undone!"):
                backup_path = Path(__file__).resolve().parent.parent / "data" / "mutation_rates_backup.csv"
                default_df = pd.read_csv(backup_path, index_col=0, header=0, names=["Marker", "Mutation rate"])
                default_df.to_csv(mutation_rates_path, index=True, header=True)
                st.success("✅ Reset to default!")
                st.rerun()

    st.markdown("")

    with st.container():
        st.markdown("##### ➕ Add New Marker")

        col_name, col_rate = st.columns([2, 1])

        with col_name:
            new_marker_name = st.text_input(
                "Marker name",
                placeholder="e.g., DYS389I",
                help="Enter a unique marker name"
            )

        with col_rate:
            new_marker_rate = st.number_input(
                "Mutation rate",
                min_value=0.0000,
                max_value=1.0000,
                value=0.0001,
                step=0.0001,
                format="%.8f",
                help="Mutation rate per generation (0-1)"
            )

        if st.button("➕ Add Marker", type="primary", width='stretch'):
            if new_marker_name and new_marker_rate is not None:
                if new_marker_name in edited_df.index:
                    st.error("❌ Marker already exists. Please choose a different name.")
                elif new_marker_rate < 0.0000 or new_marker_rate > 1.0000:
                    st.error("❌ Mutation rate must be between 0 and 1.")
                else:
                    new_row = pd.DataFrame([[new_marker_name, new_marker_rate]], columns=["Marker", "Mutation rate"])
                    new_row.set_index("Marker", inplace=True)
                    edited_df = pd.concat([edited_df, new_row], ignore_index=False)
                    edited_df.to_csv(mutation_rates_path, index=True, header=True)
                    st.success(f"✅ Marker '{new_marker_name}' added!")
                    st.rerun()
            else:
                st.error("❌ Please enter both marker name and mutation rate.")

# ==================== RIGHT COLUMN: MARKER SETS ====================
with col2:
    st.subheader("📦 Marker Sets (Kits)")

    with st.container():
        st.markdown("##### Select or Create Marker Set")

        kit_names = list(kits.keys())
        kit_names.append("➕ Create New Kit")

        selected_kit = st.selectbox(
            "Select a kit",
            kit_names,
            help="Choose an existing kit or create a new one"
        )

        if selected_kit == "➕ Create New Kit":
            st.info("💡 Enter a name for the new kit and select markers to include.")
            new_kit_name = st.text_input("Enter new kit name", placeholder="e.g., MyCustomKit")
            selected_kit = new_kit_name if new_kit_name else ""
            if selected_kit:
                kits[selected_kit] = []

        if selected_kit and selected_kit != "➕ Create New Kit":
            selected_kit_markers = kits.get(selected_kit, [])

            st.markdown("##### Configure Markers")

            marker_multiselect = st.multiselect(
                "Select markers for this kit",
                options=mutation_rates_df.index.tolist(),
                default=selected_kit_markers,
                help="Select the markers to include in the kit. Markers must exist in the mutation rates table.",
            )

            # Calculate overall mutation rate
            if marker_multiselect:
                overall_mutation_rate = 1 - math.prod([1 - mut_rate for mut_rate in mutation_rates_df.loc[marker_multiselect, "Mutation rate"]])
                st.success(f"📊 **Overall mutation rate:** {overall_mutation_rate:.6f} mutation(s) per generation\n\n({len(marker_multiselect)} markers selected)")
            else:
                st.warning("⚠️ No markers selected for this kit.")

            st.markdown("##### Actions")

            col_save_kit, col_delete_kit = st.columns([1, 1])

            with col_save_kit:
                if st.button("💾 Save Kit", type="primary", width='stretch'):
                    if selected_kit and marker_multiselect and len(marker_multiselect) > 0:
                        kits[selected_kit] = marker_multiselect
                        with open(kits_path, "w") as file:
                            json.dump(kits, file, indent=4)
                        st.success(f"✅ Kit '{selected_kit}' saved!")
                        st.rerun()
                    else:
                        st.error("❌ Please select a kit name and at least one marker.")

            with col_delete_kit:
                if st.button("🗑️ Delete Kit", type="secondary", width='stretch'):
                    if selected_kit and selected_kit in kits:
                        del kits[selected_kit]
                        with open(kits_path, "w") as file:
                            json.dump(kits, file, indent=4)
                        st.success(f"✅ Kit '{selected_kit}' deleted!")
                        st.rerun()
                    else:
                        st.error("❌ Please select a valid kit to delete.")

            st.markdown("")

            if st.button("🔄 Reset All Kits to Default", type="secondary", width='stretch', help="⚠️ This action cannot be undone!"):
                backup_path = Path(__file__).resolve().parent.parent / "data" / "kits_backup.json"
                with open(backup_path, "r") as file:
                    default_kits = json.load(file)
                with open(kits_path, "w") as file:
                    json.dump(default_kits, file, indent=4)
                st.success("✅ All kits reset to default!")
                st.rerun()