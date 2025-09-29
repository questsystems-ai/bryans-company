# app/pages/02_benchmarks.py
import re
import math
import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(page_title="5-HT2A/5-HT2B Benchmarks", layout="wide")

# ---------- Paths (adjust if your repo layout differs) ----------
P_BEST = "data/results/best_poses.parquet"               # from scripts/postprocess.py
P_CNN  = "data/results/best_cnn_by_receptor.parquet"     # from scripts/postprocess.py
P_BENCH= "benchmark_Ki_values_extended.csv"              # put next to your Streamlit entrypoint or pass full path

# ---------- Utilities ----------
def extract_numbers(s):
    """Return all numeric values (floats) found in a string (handles en/em dashes)."""
    if s is None or not isinstance(s, str) or not s.strip():
        return []
    s = s.replace("–", "-").replace("—", "-").replace("µ", "u")
    nums = re.findall(r"\d+\.?\d*", s)
    return [float(x) for x in nums]

def geom_mean(xs):
    xs = [x for x in xs if x > 0]
    if not xs:
        return None
    return float(np.exp(np.mean(np.log(xs))))

def parse_ki_nm(ki_str):
    """
    Parse Ki text into a single representative Ki in nM.
    Strategy: collect all numbers, take geometric mean; if only one, return it.
    """
    xs = extract_numbers(ki_str)
    if not xs:
        return None
    return geom_mean(xs)

def nm_to_pKi(ki_nm):
    """Convert nM -> pKi (Molar): pKi = -log10(Ki[M]); Ki[nM] -> Ki[M] = Ki*1e-9"""
    if ki_nm is None or ki_nm <= 0:
        return None
    return -math.log10(ki_nm * 1e-9)

# ---------- Load data ----------
@st.cache_data
def load_gnina():
    try:
        best = pd.read_parquet(P_BEST)
    except Exception:
        best = pd.DataFrame()
    try:
        piv = pd.read_parquet(P_CNN)  # index=ligand, columns=receptor, values=cnn
    except Exception:
        piv = pd.DataFrame()
    return best, piv

@st.cache_data
def load_benchmark(path):
    df = pd.read_csv(path)
    # Normalize ligand names to a canonical form for joining
    df["ligand_key"] = df["Ligand"].str.strip().str.lower()
    # parse Ki strings -> numeric (nM)
    df["Ki_5HT2A_nM_num"] = df["Ki_5HT2A_nM"].apply(parse_ki_nm)
    df["Ki_5HT2B_nM_num"] = df["Ki_5HT2B_nM"].apply(parse_ki_nm)
    # compute pKi
    df["pKi_5HT2A"] = df["Ki_5HT2A_nM_num"].apply(nm_to_pKi)
    df["pKi_5HT2B"] = df["Ki_5HT2B_nM_num"].apply(nm_to_pKi)
    return df

best, piv = load_gnina()
bench = load_benchmark(P_BENCH)

st.title("Benchmarking Gnina vs Experimental Ki (5-HT2A / 5-HT2B)")

# Info panel
with st.expander("What this page does", expanded=False):
    st.markdown(
        "- Loads **Gnina CNN scores** per ligand & receptor from your pipeline.\n"
        "- Loads a small **benchmark Ki set** (psilocin, DMT, LSD, serotonin).\n"
        "- Converts Ki ranges → numeric estimates (geometric mean), computes **pKi**.\n"
        "- Shows **tables + scatter plots** for 5-HT2A, 5-HT2B, and selectivity (ΔpKi vs ΔCNN)."
    )

# Early checks
if piv.empty:
    st.error(f"Could not find {P_CNN}. Run your docking + postprocess first.")
    st.stop()

# Make ligand key in Gnina matrix
piv = piv.copy()
piv.index.name = "ligand"
piv.reset_index(inplace=True)
piv["ligand_key"] = piv["ligand"].str.strip().str.lower()

# Heuristic mapping of receptor column names
# Accept common labels: '6WGT','7WC5','4IB4','5HT2A','5HT2B', etc.
def pick_col(candidates, cols):
    for c in candidates:
        if c in cols: return c
    # fuzzy contains
    for c in cols:
        if any(tok.lower() in c.lower() for tok in candidates):
            return c
    return None

cols = list(piv.columns)
col_a = pick_col(["6WGT", "7WC5", "5HT2A", "2a", "ht2a"], cols)
col_b = pick_col(["4IB4", "5HT2B", "2b", "ht2b"], cols)

if not col_a or not col_b:
    st.warning(
        f"Could not auto-detect 5-HT2A/5-HT2B columns from {cols}. "
        "Use column names that include '6WGT'/'7WC5' or '5HT2A' for A, and '4IB4' or '5HT2B' for B."
    )

# Merge benchmarks with cnn per ligand
merged = pd.merge(bench, piv[["ligand_key", "ligand", col_a, col_b]],
                  on="ligand_key", how="left", suffixes=("", "_cnn"))

# Derived columns
merged = merged.rename(columns={col_a: "CNN_5HT2A", col_b: "CNN_5HT2B"})
merged["ΔCNN (A-B)"] = merged["CNN_5HT2A"] - merged["CNN_5HT2B"]
merged["ΔpKi (A-B)"] = merged["pKi_5HT2A"] - merged["pKi_5HT2B"]

# Display
left, right = st.columns([1.1, 1])
with left:
    st.subheader("Benchmark table (parsed Ki → pKi) joined with your CNN scores")
    st.dataframe(
        merged[[
            "Ligand", "Notes",
            "Ki_5HT2A_nM", "Ki_5HT2A_nM_num", "pKi_5HT2A",
            "Ki_5HT2B_nM", "Ki_5HT2B_nM_num", "pKi_5HT2B",
            "CNN_5HT2A", "CNN_5HT2B", "ΔpKi (A-B)", "ΔCNN (A-B)"
        ]].sort_values("Ligand"),
        use_container_width=True, height=420
    )

with right:
    st.markdown("**Legend**")
    st.markdown(
        "- **pKi** = −log10(Ki in M). Higher pKi = tighter binding.\n"
        "- **CNN** = Gnina CNNscore (higher = better predicted pose/affinity).\n"
        "- **Δ** = A − B. Positive ΔpKi means more selective for **5-HT2A** over **5-HT2B**."
    )

# Plots
st.markdown("---")
st.subheader("Correlation plots")

import altair as alt

def scatter(df, x, y, color="Ligand", tooltip=None, title=""):
    tooltip = tooltip or [alt.Tooltip(c) for c in df.columns]
    chart = alt.Chart(df).mark_circle(size=120).encode(
        x=alt.X(x, scale=alt.Scale(zero=False)),
        y=alt.Y(y, scale=alt.Scale(zero=False)),
        color=alt.Color(color, legend=alt.Legend(title="Ligand")),
        tooltip=tooltip
    ).properties(height=340, title=title).interactive()
    return chart

valid_a = merged.dropna(subset=["pKi_5HT2A", "CNN_5HT2A"])
valid_b = merged.dropna(subset=["pKi_5HT2B", "CNN_5HT2B"])
valid_sel = merged.dropna(subset=["ΔpKi (A-B)", "ΔCNN (A-B)"])

c1 = scatter(
    valid_a, "pKi_5HT2A:Q", "CNN_5HT2A:Q",
    title="5-HT2A: CNNscore vs pKi (higher is better on both axes)"
)
c2 = scatter(
    valid_b, "pKi_5HT2B:Q", "CNN_5HT2B:Q",
    title="5-HT2B: CNNscore vs pKi (ideally lower binding than 2A)"
)
c3 = scatter(
    valid_sel, "ΔpKi (A-B):Q", "ΔCNN (A-B):Q",
    title="Selectivity: ΔCNN vs ΔpKi (A−B)  (upper-right is best for 2A-selective)"
)

s1, s2 = st.columns(2)
with s1: st.altair_chart(c1, use_container_width=True)
with s2: st.altair_chart(c2, use_container_width=True)

st.altair_chart(c3, use_container_width=True)

# Download merged table
csv_out = merged.to_csv(index=False).encode()
st.download_button(
    "⬇️ Download merged benchmark + CNN table (CSV)",
    csv_out,
    file_name="benchmarks_vs_cnn.csv",
    mime="text/csv"
)

st.caption(
    "Notes:\n"
    "- Ki values are literature aggregates; some entries reflect ranges and species differences. "
    "We parse a representative **geometric mean** for plotting.\n"
    "- For more precise work, replace the CSV with exact human 5-HT2A/2B Ki entries from PDSP/ChEMBL/BindingDB."
)
