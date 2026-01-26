import streamlit as st
from pedigree_lr.visualization import st_visualize_pedigree
from pedigree_lr.models import Pedigree
from configparser import ConfigParser
from pathlib import Path

st.set_page_config(
    page_title="Pedigree Builder",
    page_icon="🏗️",
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

st.title("🏗️ Pedigree Builder")
st.markdown("Build Y-chromosome pedigrees interactively by adding individuals and defining relationships.")
st.markdown("---")

# Initialize config
config_path = Path(__file__).resolve().parent.parent / "data" / "config.ini"
if "global_config" not in st.session_state:
    global_config = ConfigParser()
    global_config.optionxform = str  # type: ignore
    global_config.read(config_path)
    st.session_state.global_config = global_config

# Initialize pedigree
if 'pedigree' not in st.session_state or st.session_state.pedigree is None:
    st.session_state["pedigree"] = Pedigree()

individuals = [i for i in st.session_state.pedigree.individuals]

with st.sidebar:
    st.markdown("### 👤 Add Individual")

    if len(individuals) > 0:
        st.markdown("##### Select Father (optional)")
        father = st.selectbox(
            "Father",
            list(individuals),
            index=len(individuals) - 1,
            format_func=lambda x: x.name,
            help="Select the father for the new individual. Leave as default for founders.",
            label_visibility="collapsed"
        )
    else:
        father = None
        st.info("💡 Start by adding the founder (oldest generation)")

    new_individual = st.text_input(
        "Individual name/ID",
        value="",
        help="Enter a unique name or ID for the individual.",
        placeholder="e.g., John_Doe or ID001",
        key="new_individual_key"
    ).strip()

    if st.button("➕ Add Individual", type="primary", width='stretch'):
        if new_individual and new_individual != "":
            if new_individual in [i.name for i in individuals]:
                st.error(f"❌ Individual '{new_individual}' already exists. Use a unique name/ID.")
            else:
                st.session_state.pedigree.add_individual(
                    individual_id=new_individual,
                    name=new_individual
                )
                if father and father.id != "None":
                    st.session_state.pedigree.add_relationship(
                        parent_id=father.id,
                        child_id=new_individual
                    )
                st.success(f"✅ Added: {new_individual}")
                st.rerun()
        else:
            st.error("❌ Please enter a name/ID for the individual.")

    if len(individuals) > 0:
        st.divider()
        st.markdown("### 🗑️ Remove Individual")

        individual_to_remove = st.selectbox(
            "Select individual to remove",
            list(individuals),
            format_func=lambda x: x.name,
            help="Removes the selected individual and all descendants."
        )

        if st.button("🗑️ Remove Individual", type="secondary", width='stretch'):
            st.session_state.pedigree.remove_individual(individual_to_remove.id)
            st.success(f"✅ Removed: {individual_to_remove.name}")
            st.rerun()

        st.divider()
        st.markdown("### 💾 Save Pedigree")

        st.info("💡 Download the pedigree file to use it in simulations.")
        st.download_button(
            "📥 Download Pedigree (TGF)",
            st.session_state.pedigree.write_to_tgf(),
            file_name="pedigree.tgf",
            help="Download pedigree in TGF format",
            width='stretch'
        )

# Main content area
if len(individuals) > 0:
    with st.container():
        st.subheader("📊 Pedigree Visualization")
        st.markdown(f"**Total individuals:** {len(individuals)}")

        st_visualize_pedigree(
            st.session_state.pedigree,
            global_config=st.session_state.global_config
        )
else:
    st.session_state.pedigree = None
    with st.container():
        st.info("### 👋 Welcome to the Pedigree Builder!")
        st.markdown("""
        This tool allows you to build Y-chromosome pedigrees interactively.

        **How to use:**
        1. 📝 **Add the founder** - Start with the oldest generation (no father selected)
        2. 👨‍👦 **Add children** - Select a father and add descendants
        3. 🔄 **Repeat** - Continue building the pedigree generation by generation
        4. 💾 **Download** - Save your pedigree as a TGF file

        **Tips:**
        - Build from oldest to youngest generation
        - Use unique names or IDs for each individual
        - Removing an individual also removes all descendants
        - Download the pedigree to use in the Home page
        """)

        col_tip1, col_tip2, col_tip3 = st.columns(3)
        with col_tip1:
            st.markdown("##### 🎯 Start Simple")
            st.markdown("Begin with 2-3 generations")
        with col_tip2:
            st.markdown("##### 📏 Y-Chromosome")
            st.markdown("Only male lineage is tracked")
        with col_tip3:
            st.markdown("##### 📝 Naming")
            st.markdown("Use clear, unique identifiers")