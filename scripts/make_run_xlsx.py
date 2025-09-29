#!/usr/bin/env python3
import os, json, subprocess, sys, math, tempfile
from pathlib import Path
from datetime import datetime

import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import Draw

import xlsxwriter

REPO = Path.cwd()
DEFAULT_TAG = "try1"
TARGETS = ["5HT2A_active_7WC5", "5HT2A_active_6WGT", "5HT2B_4IB4"]

# -----------------------
# helpers
# -----------------------
def nvidia_info():
    try:
        out = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=3)
        if out.returncode == 0:
            line = next((l for l in out.stdout.splitlines() if "NVIDIA-SMI" in l and "Driver Version" in l), "")
            return line.strip()
    except Exception:
        pass
    return "nvidia-smi unavailable"

def rd_img_from_mol(m, size=(320,240)):
    try:
        return Draw.MolToImage(m, size=size)
    except Exception:
        return None

def load_tables(tag):
    base = REPO / f"data/screens/{tag}"
    raw = pd.read_csv(base/"raw_best_by_target.csv")   # ligand,target,affinity,cnnpose,cnnaff,log
    # primary summary (best_2A, best_2B, delta) + CNN if available
    try:
        summ = pd.read_csv(base/"summary_with_cnn.csv", index_col=0)
    except Exception:
        summ = pd.read_csv(base/"summary_per_ligand.csv", index_col=0)
        # add cnn_2A / cnn_2B from raw
        cnn_2a = (raw[raw.target.str.contains("5HT2A_active_")]
                  .sort_values(["ligand","cnnpose"], ascending=[True,False])
                  .drop_duplicates("ligand")[["ligand","cnnpose"]]
                  .rename(columns={"cnnpose":"cnn_2A"}).set_index("ligand"))
        cnn_2b = (raw[raw.target=="5HT2B_4IB4"]
                  .sort_values(["ligand","cnnpose"], ascending=[True,False])
                  .drop_duplicates("ligand")[["ligand","cnnpose"]]
                  .rename(columns={"cnnpose":"cnn_2B"}).set_index("ligand"))
        summ = summ.join(cnn_2a, how="left").join(cnn_2b, how="left")
    return raw, summ

def load_ligands(tag):
    # Prefer one-SDF-per-ligand under screen tag
    lig_dir = REPO / f"data/screens/{tag}/ligands"
    ligs={}
    if lig_dir.exists():
        for sdf in lig_dir.glob("*.sdf"):
            sup = Chem.SDMolSupplier(str(sdf), removeHs=False)
            m = next((mol for mol in sup if mol), None)
            if m:
                nm = sdf.stem
                m.SetProp("_Name", nm)
                ligs[nm] = m
        if ligs:
            return ligs
    # Fallback: multi SDF with _Name
    multi = REPO / "data/ligands/lib.sdf"
    if multi.exists():
        for m in Chem.SDMolSupplier(str(multi), removeHs=False):
            if not m: continue
            nm = m.GetProp("_Name") if m.HasProp("_Name") else f"lig_{len(ligs)+1}"
            m.SetProp("_Name", nm)
            ligs[nm] = m
    return ligs

def smiles_for_mol(m):
    try:
        return Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True)
    except Exception:
        return ""

def pretty_name(lig_name: str, smi: str | None = None) -> str:
            overrides = {
                "Serotonin": "Serotonin",
                "LSD": "Lysergic acid diethylamide",
                "Psilocin": "Psilocin",
                "DMT_like": "N,N-Dimethyltryptamine",
                "MET_like": "N-Methyl-N-ethyltryptamine",
                "TRP_like": "Tryptamine",
            }
            if lig_name in overrides: return overrides[lig_name]
            low = lig_name.lower()
            if "dmt" in low: return "N,N-Dimethyltryptamine"
            if "met" in low: return "N-Methyl-N-ethyltryptamine"
            if "trp" in low or "trypt" in low: return "Tryptamine"
            return lig_name

# -----------------------
# writer
# -----------------------
def write_xlsx(tag, raw, summ, lig_mols, grids, params):
    out_xlsx = REPO / f"data/screens/{tag}/run_report_{tag}.xlsx"
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb = xlsxwriter.Workbook(str(out_xlsx))
    img_tmp = Path(tempfile.gettempdir())

    # formats
    fmt_hdr = wb.add_format({"bold": True, "bg_color": "#F2F2F2", "border":1})
    fmt_cell = wb.add_format({"border":1})
    fmt_num  = wb.add_format({"border":1, "num_format":"0.00"})
    fmt_small = wb.add_format({"font_size":9})
    fmt_title = wb.add_format({"bold": True, "font_size": 14})

    # --- Params sheet
    wsP = wb.add_worksheet("Params")
    wsP.write(0,0,f"5-HT Docking Run Report — {tag}", fmt_title)
    wsP.write(1,0,"Generated")
    wsP.write(1,1,datetime.now().strftime('%Y-%m-%d %H:%M'))
    wsP.write(2,0,"GPU/driver")
    wsP.write(2,1,nvidia_info())
    wsP.write(4,0,"Exhaustiveness"); wsP.write_number(4,1, params["exh"])
    wsP.write(5,0,"num_modes"); wsP.write_number(5,1, params["modes"])
    wsP.write(6,0,"Box size (Å)"); wsP.write_number(6,1, params["size"])
    wsP.write(7,0,"Protonation"); wsP.write(7,1,"Dimorphite-DL pH 6.8–8.0 (RDKit 3D)")
    row=9
    for tgt,meta in grids.items():
        wsP.write(row,0,f"Grid center — {tgt}")
        wsP.write(row,1,", ".join([f"{x:.2f}" for x in meta["center"]]))
        row+=1
    wsP.set_column(0,0,28); wsP.set_column(1,1,60)

    # --- Summary sheet
    ws = wb.add_worksheet("Summary")
    headers = ["Ligand", "Friendly Name", "Structure", "SMILES",
               "best_2A (kcal/mol)", "best_2B (kcal/mol)", "Δ(2A–2B)", "cnn_2A", "cnn_2B"]
    for j,h in enumerate(headers): ws.write(0,j,h,fmt_hdr)
    ws.set_row(0, 22)
    ws.set_column(0,0,18)  # Ligand
    ws.set_column(1,1,26)  # Friendly Name
    ws.set_column(2,2,24)  # Structure
    ws.set_column(3,3,42)  # SMILES
    ws.set_column(4,8,16)  # numeric cols

    # order ligands
    view = summ.copy()
    view["rank_key"] = list(zip(view["delta_2A_2B"], -view["best_2A"]))
    view = view.sort_values(["delta_2A_2B","best_2A"])

    r=1
    for lig, rowV in view.iterrows():
        # ligand name & smiles
        m = lig_mols.get(lig, None)
        smi = smiles_for_mol(m) if m else ""
        friendly = pretty_name(lig, smi)

        ws.write(r,0, lig, fmt_cell)
        ws.write(r,1, friendly, fmt_cell)
        ws.write(r,3, smi, fmt_cell)

        # numbers
        for col_key, col_idx in [("best_2A",4), ("best_2B",5), ("delta_2A_2B",6),
                                 ("cnn_2A",7), ("cnn_2B",8)]:
            val = rowV.get(col_key, None)
            if pd.notna(val):
                ws.write_number(r, col_idx, float(val), fmt_num)
            else:
                ws.write(r, col_idx, "NA", fmt_cell)

        # image in cell
        if m:
            img = rd_img_from_mol(m, size=(320,240))
            if img:
                tmp_png = img_tmp / f"{lig}.png"
                try:
                    img.save(tmp_png)
                    # row height for image; adjust scaling slightly
                    ws.set_row(r, 190)
                    ws.insert_image(r, 2, str(tmp_png),
                                    {"x_scale":0.9, "y_scale":0.9, "object_position":1})
                except Exception:
                    pass

        r += 1

    # --- ByTarget sheet
    ws2 = wb.add_worksheet("ByTarget")
    raw = pd.read_csv(REPO/f"data/screens/{DEFAULT_TAG}/raw_best_by_target.csv")
    cols = ["ligand","target","affinity","cnnpose","cnnaff","log"]
    for j,h in enumerate(cols): ws2.write(0,j,h,fmt_hdr)
    ws2.set_row(0,22); ws2.set_column(0,0,22); ws2.set_column(1,1,22)
    ws2.set_column(2,4,16); ws2.set_column(5,5,72)
    rr=1
    for t in raw.itertuples():
        ws2.write(rr,0,t.ligand,fmt_cell)
        ws2.write(rr,1,t.target,fmt_cell)
        try: ws2.write_number(rr,2,float(t.affinity),fmt_num)
        except: ws2.write(rr,2,"NA",fmt_cell)
        try: ws2.write_number(rr,3,float(t.cnnpose),fmt_num)
        except: ws2.write(rr,3,"NA",fmt_cell)
        try: ws2.write_number(rr,4,float(t.cnnaff),fmt_num)
        except: ws2.write(rr,4,"NA",fmt_cell)
        ws2.write(rr,5,t.log,fmt_cell)
        rr+=1

    wb.close()
    return out_xlsx

# -----------------------
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default=DEFAULT_TAG)
    ap.add_argument("--size", type=int, default=24)
    ap.add_argument("--exh", type=int, default=128)
    ap.add_argument("--modes", type=int, default=30)
    args = ap.parse_args()

    raw, summ = load_tables(args.tag)
    grids = json.loads((REPO/"data/grids/grids.json").read_text())
    lig_mols = load_ligands(args.tag)
    out = write_xlsx(args.tag, raw, summ, lig_mols, grids, vars(args))
    print(f"[done] wrote {out}")

if __name__ == "__main__":
    main()
