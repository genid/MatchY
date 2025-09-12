import streamlit as st
from configparser import ConfigParser
from pathlib import Path

st.set_page_config(
    page_title="Software settings",
    page_icon="⚙️",
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

st.sidebar.header("Software settings")

config_path = Path(__file__).resolve().parent.parent / "data" / "config.ini"
global_config = st.session_state.global_config

st.info("Edit the software settings. The changes will be applied to the current session.")

for section in global_config.sections():
    st.subheader(section)
    for key, value in global_config.items(section):
        new_value = st.text_input(f"{key}:", value=value, key=key)
        if new_value != value:
            global_config.set(section, key, str(new_value))

st.info("Click 'Save changes' to apply the changes to the software settings for future sessions.")
if st.button("Save changes", type="primary"):
    with open(config_path, "w") as configfile:
        global_config.write(configfile)
    st.success("Changes saved successfully.")
    st.session_state.global_config = global_config