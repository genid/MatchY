import streamlit as st
from configparser import ConfigParser
from pathlib import Path

st.set_page_config(
    page_title="Software settings",
    page_icon="⚙️",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.logo("assets/logo_minimal.png", icon_image="assets/icon.png")
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

st.title("⚙️ Software Settings")
st.markdown("Configure application settings and parameters for simulations.")
st.markdown("---")

config_path = Path(__file__).resolve().parent.parent / "data" / "config.ini"

# Initialize global_config if it doesn't exist in session state
if "global_config" not in st.session_state:
    global_config = ConfigParser()
    global_config.optionxform = str  # type: ignore
    global_config.read(config_path)
    st.session_state.global_config = global_config
else:
    global_config = st.session_state.global_config

st.info("💡 Edit the settings below. Changes will be applied to the current session immediately and can be saved for future sessions.")

# Group settings by section
for section in global_config.sections():
    with st.container():
        st.subheader(f"📋 {section}")

        # Create columns for better layout
        num_items = len(list(global_config.items(section)))
        cols_per_row = 2

        items = list(global_config.items(section))
        for i in range(0, num_items, cols_per_row):
            cols = st.columns(cols_per_row)
            for j, (key, value) in enumerate(items[i:i+cols_per_row]):
                with cols[j]:
                    new_value = st.text_input(
                        f"{key}",
                        value=value,
                        key=key,
                        help=f"Current value: {value}"
                    )
                    if new_value != value:
                        global_config.set(section, key, str(new_value))

        st.markdown("")  # Spacing

st.markdown("---")

# Save section
with st.container():
    st.subheader("💾 Save Settings")
    st.info("💡 Click the button below to save your changes to the configuration file for future sessions.")

    col_save, col_spacer = st.columns([1, 3])

    with col_save:
        if st.button("💾 Save Changes", type="primary", width='stretch'):
            with open(config_path, "w") as configfile:
                global_config.write(configfile)
            st.session_state.global_config = global_config
            st.success("✅ Settings saved successfully!")