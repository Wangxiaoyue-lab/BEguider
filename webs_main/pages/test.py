import streamlit as st

# Create a state variable to keep track of the currently open expander
state = st.session_state
if "open_expander" not in state:
    state.open_expander = None

# Create multiple expanders
with st.expander("Expander 1", expanded=state.open_expander == "expander1"):
    st.write("Content of expander 1")
    if st.button("Open next expander"):
        state.open_expander = "expander2"

with st.expander("Expander 2", expanded=state.open_expander == "expander2"):
    st.write("Content of expander 2")
    if st.button("Open next expander"):
        state.open_expander = "expander3"

with st.expander("Expander 3", expanded=state.open_expander == "expander3"):
    st.write("Content of expander 3")
