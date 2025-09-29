import streamlit as st
import pandas as pd
import py3Dmol
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
st.set_page_config(page_title="5-HT Docking Dashboard", layout="wide")
st.title("5-HT Docking Dashboard (GNINA + CNN rescoring)")

best_path = BASE/"data"/"results"/"best_poses.parquet"
if not best_path.exists():
    st.warning("No results yet. Run docking and postprocess first.")
    st.stop()
df = pd.read_parquet(best_path)

ligands = sorted(df["ligand"].unique())
receptors = sorted(df["receptor"].unique())

col1, col2, col3 = st.columns([1,1,2])
with col1:
    lig = st.selectbox("Ligand", ligands)
with col2:
    rec = st.selectbox("Receptor", receptors)
subset = df[(df.ligand==lig) & (df.receptor==rec)]
if subset.empty:
    st.warning("No pose for this selection.")
    st.stop()
vina = subset["vina"].values[0]
cnn = subset["cnn"].values[0]
st.markdown(f"**CNNscore:** {cnn:.3f} &nbsp;&nbsp; **Vina:** {vina:.2f}")

sdf_path = BASE/"data"/"results"/f"dock_{rec}"/"poses.sdf"
viewer = py3Dmol.view(width=900, height=600)
viewer.addModel(open(sdf_path).read(), "sdf")
viewer.setStyle({"stick": {}})
viewer.zoomTo()
st.components.v1.html(viewer._make_html(), height=620, scrolling=False)

st.caption("Sort/filter below and export as needed.")
st.dataframe(df.sort_values(["receptor","cnn"], ascending=[True,False]))
