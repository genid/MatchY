import streamlit as st
from pedigree_lr.visualization import st_visualize_pedigree
from pedigree_lr.models import Pedigree
from hashlib import sha256

st.set_page_config(
    page_title="Pedigree Builder",
    page_icon="🏗️",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.sidebar.header("Pedigree Builder")

if 'pedigree_builder' not in st.session_state:
    st.session_state["pedigree_builder"] = Pedigree()

individuals = [i.id for i in getattr(st.session_state.pedigree_builder, 'individuals', [])]

with st.sidebar:
    if len(individuals) > 0:
        father = st.selectbox("Select the father:", list(individuals), index=len(individuals) - 1)
    else:
        father = None

    new_individual = st.text_input("Enter the unique name/ID of the individual:",
                                   value="",
                                   help="Start with the oldest generation.",
                                   key="new_individual_key"
                                   ).replace(" ", "_")

    if st.button("Add Individual"):
        if new_individual:
            if new_individual in individuals:
                st.error(f"Individual {new_individual} already exists. Use a unique name/ID.")
            else:
                st.session_state.pedigree_builder.add_individual(new_individual, new_individual)
                if father and father != "None":
                    st.session_state.pedigree_builder.add_relationship(father, new_individual)
                st.success(f"Added individual: {new_individual}")
                st.session_state.update({"add_individual_triggered": False})
                st.rerun()
        else:
            st.error("Please enter a name/ID for the individual.")

    # add option to remove individuals (and all his descendants) from the pedigree
    if len(individuals) > 0:
        st.markdown("---")
        individual_to_remove = st.selectbox("Select the individual to remove:", list(individuals))
        if st.button("Remove Individual"):
            st.session_state.pedigree_builder.remove_individual(individual_to_remove)
            st.success(f"Removed individual: {individual_to_remove}.")
            st.rerun()

        st.markdown("---")
        st.download_button("Download Pedigree", st.session_state.pedigree_builder.write_to_tgf(), file_name="pedigree.tgf")
        sha_hash = sha256(st.session_state.pedigree_builder.write_to_tgf()).hexdigest()
        st.write(f"SHA-256 hash: {sha_hash}")

if len(individuals) > 0:
    st.write("Pedigree:")
    st_visualize_pedigree(st.session_state.pedigree_builder)
