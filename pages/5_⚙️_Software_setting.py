import streamlit as st
from configparser import ConfigParser
from pathlib import Path

st.set_page_config(
    page_title="Software settings",
    page_icon="⚙️",
    layout="wide",
    initial_sidebar_state="expanded",
)

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