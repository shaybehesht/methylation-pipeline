import streamlit as st

from app.state import initialize
from core.thresholds import REGISTRY, validate, widget_values

initialize()
st.title("🥭 3. Thresholds")
st.write("Every control and its explanation come from the same registry.")

values = {}
for key, spec in REGISTRY.items():
    current = st.session_state.thresholds.get(key, spec.default)
    minimum, maximum, current, step = widget_values(spec, current)
    if spec.widget == "slider":
        values[key] = st.slider(
            spec.label, minimum, maximum, current, step, key=f"threshold_{key}"
        )
    else:
        values[key] = st.number_input(
            spec.label, min_value=minimum, max_value=maximum,
            value=current, step=step, key=f"threshold_{key}",
        )
    with st.expander("Why this default?"):
        st.write(spec.rationale)
        st.caption(f"Default: {spec.default}")

st.session_state.thresholds = validate(values)
