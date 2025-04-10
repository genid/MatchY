import streamlit as st
from pedigree_lr.visualization import st_visualize_pedigree
from pedigree_lr.models import Pedigree

st.set_page_config(
    page_title="Pedigree Builder",
    page_icon="🏗️",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.sidebar.header("Pedigree Builder")

if 'pedigree' in st.session_state and st.session_state.pedigree is not None:
    st.info("Since you have uploaded a pedigree, you can change it here.")
    st.session_state.pedigree_builder = st.session_state.pedigree
    st.session_state.pedigree = None

if 'pedigree_builder' not in st.session_state:
    st.session_state["pedigree_builder"] = Pedigree()

individuals = [i for i in getattr(st.session_state.pedigree_builder, 'individuals', [])]

with st.sidebar:
    if len(individuals) > 0:
        father = st.selectbox("Select the father:",
                              list(individuals),
                              index=len(individuals) - 1,
                              format_func=lambda x: x.name
                              )
    else:
        father = None

    st.info("Add individuals to the pedigree. Start with the oldest generation.")
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
                st.session_state.pedigree_builder.add_individual(individual_id=new_individual,
                                                                 name=new_individual)
                if father and father.id != "None":
                    st.session_state.pedigree_builder.add_relationship(parent_id=father.id,
                                                                       child_id=new_individual)
                st.success(f"Added individual: {new_individual}")
                st.session_state.update({"add_individual_triggered": False})
                st.rerun()
        else:
            st.error("Please enter a name/ID for the individual.")

    # TODO: add option to remove individuals (and all his descendants) from the pedigree
    if len(individuals) > 0:
        st.divider()
        individual_to_remove = st.selectbox("Select the individual to remove:",
                                            list(individuals),
                                            format_func=lambda x: x.name)
        if st.button("Remove Individual"):
            st.session_state.pedigree_builder.remove_individual(individual_to_remove.id)
            st.success(f"Removed individual: {individual_to_remove.name}.")
            st.rerun()

        st.download_button("Download Pedigree", st.session_state.pedigree_builder.write_to_tgf(), file_name="pedigree.tgf", type="tertiary")

        if st.button("Load directly to simulation",
                     type="primary"):
            st.session_state.pedigree = st.session_state.pedigree_builder
            st.success("Pedigree loaded to simulation.")
            st.info("Go to Home to start the simulation.")
            st.warning("Remember to download the pedigree file if you want to use it for future simulations.")

if len(individuals) > 0:
    st_visualize_pedigree(st.session_state.pedigree_builder,
                          global_config=st.session_state.global_config)
