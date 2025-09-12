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

st.sidebar.header("Pedigree Builder")

config_path = Path(__file__).resolve().parent.parent / "data" / "config.ini"
global_config = ConfigParser()
global_config.optionxform = str  # type: ignore
global_config.read(config_path)
st.session_state.global_config = global_config

if 'pedigree' not in st.session_state or st.session_state.pedigree is None:
    st.session_state["pedigree"] = Pedigree()

individuals = [i for i in st.session_state.pedigree.individuals]

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

    if st.button("Add Individual",
                 type="primary"):
        if new_individual:
            if new_individual in [i.name for i in individuals]:
                st.error(f"Individual {new_individual} already exists. Use a unique name/ID.")
            else:
                st.session_state.pedigree.add_individual(individual_id=new_individual,
                                                                 name=new_individual)
                if father and father.id != "None":
                    st.session_state.pedigree.add_relationship(parent_id=father.id,
                                                                       child_id=new_individual)
                st.success(f"Added individual: {new_individual}")
                st.rerun()
        else:
            st.error("Please enter a name/ID for the individual.")

    if len(individuals) > 0:
        individual_to_remove = st.selectbox("Select the individual (and all his descendants) to remove:",
                                            list(individuals),
                                            format_func=lambda x: x.name)
        if st.button("Remove individual(s)",
                     type="secondary"):
            st.session_state.pedigree.remove_individual(individual_to_remove.id)
            st.success(f"Removed individual: {individual_to_remove.name}.")
            st.rerun()

        st.warning("Remember to download the pedigree file if you want to use it for future simulations.")
        st.download_button("Download Pedigree",
                           st.session_state.pedigree.write_to_tgf(),
                           file_name="pedigree.tgf",
                           type="secondary")

if len(individuals) > 0:
    st_visualize_pedigree(st.session_state.pedigree,
                          global_config=st.session_state.global_config)
else:
    st.session_state.pedigree = None