#!/usr/bin/env python3
import os, json, subprocess, sys, math
from pathlib import Path
import tempfile
from datetime import datetime
import pandas as pd
from tqdm import tqdm

# drawing
from rdkit import Chem
from rdkit.Chem import Draw

# PDF
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4, landscape
from reportlab.lib.units import cm
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle, SimpleDocTemplate, Image, Spacer
from reportlab.lib.styles import getSampleStyleSheet

HERE = Path.cwd()
REPO = HERE
DEFAULT_TAG = "try1"

def rd_img_from_mol(m, size=(220, 180)):
    """Return a PIL image for a molecule."""
    try:
        img = Draw.MolToImage(m, size=size)
        return img
    except Exception:
        return None

def nvidia_info():
    try:
        out = subprocess.run(["nvidia-smi"], capture_output=True, text=True, timeout=3)
        if out.returncode == 0:
            lines = out.stdout.splitlines()
            # first line has “NVIDIA-SMI X  Driver Version: Y  CUDA Version: Z”
            hdr = next((l for l in lines if "NVIDIA-SMI" in l and "Driver Version" in l), "")
            return hdr.strip()
    except Exception:
        pass
    return "nvidia-smi unavailable"

def load_tables(tag):
    base = REPO / f"data/screens/{tag}"
    raw = pd.read_csv(base/"raw_best_by_target.csv")
    summ = pd.read_csv(base/"summary_per_ligand.csv", index_col=0)
    # augment with cnn_2A / cnn_2B if available
    try:
        swc = pd.read_csv(base/"summary_with_cnn.csv", index_col=0)
        summ = swc
    except Exception:
        # compute quick cnn columns if missing
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
    lig_dir = REPO / f"data/screens/{tag}/ligands"
    if lig_dir.exists():
        # one SDF per ligand
        ligs = {}
        for sdf in lig_dir.glob("*.sdf"):
            suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
            m = next((mol for mol in suppl if mol is not None), None)
            if m is not None:
                nm = sdf.stem
                m.SetProp("_Name", nm)
                ligs[nm] = m
        if ligs:
            return ligs
    # fallback: multi SDF (names in _Name)
    multi = REPO / "data/ligands/lib.sdf"
    ligs={}
    if multi.exists():
        for m in Chem.SDMolSupplier(str(multi), removeHs=False):
            if m is None: continue
            nm = m.GetProp("_Name") if m.HasProp("_Name") else f"lig_{len(ligs)+1}"
            m.SetProp("_Name", nm)
            ligs[nm]=m
    return ligs

def page_header(c, title):
    c.setFont("Helvetica-Bold", 16)
    c.drawString(2*cm, 19*cm, title)
    c.setFont("Helvetica", 10)
    c.drawString(2*cm, 18.4*cm, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    c.drawString(2*cm, 18.0*cm, nvidia_info())
    c.line(2*cm, 17.8*cm, 28*cm, 17.8*cm)

def add_params_page(pdf_path, tag, grids, params):
    c = canvas.Canvas(str(pdf_path), pagesize=landscape(A4))
    page_header(c, f"5-HT Docking Run Report — {tag} (Parameters)")

    styles = getSampleStyleSheet()
    story = []

    # parameter listing
    rows = [
        ["Receptors", ", ".join(["6WGT (2A active)", "7WC5 (2A active)", "4IB4 (2B)"])],
        ["Box size (Å)", str(params["size"])],
        ["Exhaustiveness", str(params["exh"])],
        ["num_modes", str(params["modes"])],
        ["Protonation", "Dimorphite-DL, pH window 6.8–8.0 (3D by RDKit)"],
        ["CNN scoring", "Rescore mode"],
        ["GPU/driver", nvidia_info()],
    ]
    # grid centers
    for tgt, meta in grids.items():
        ctr = ", ".join([f"{x:.2f}" for x in meta["center"]])
        rows.append([f"Grid center — {tgt}", ctr])

    t = Table(rows, colWidths=[6*cm, 22*cm])
    t.setStyle(TableStyle([
        ("FONT", (0,0), (-1,-1), "Helvetica", 10),
        ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke),
        ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
    ]))
    w,h = t.wrapOn(c, 26*cm, 14*cm)
    t.drawOn(c, 2*cm, 17*cm - h)

    # literature note
    txt = ("<para fontname=Helvetica size=9>"
           "<b>Literature baseline:</b> Docking affinities for strong 5-HT2A agonists "
           "typically fall around ~–8 to –10 kcal/mol on active-state structures. "
           "Per-ligand literature docking values are not provided here and can be added later."
           "</para>")
    p = Paragraph(txt, styles["Normal"])
    w,h = p.wrap(26*cm, 4*cm)
    p.drawOn(c, 2*cm, 2.5*cm)

    c.showPage()
    c.save()

def add_table_pages_append(pdf_path, df_view, lig_mols, params):
    # Reopen and append pages
    c = canvas.Canvas(str(pdf_path), pagesize=landscape(A4))
    c.setAuthor("5ht-platform-starter")
    styles = getSampleStyleSheet()

    # chunk ligands per page
    rows_per_page = 6
    items = list(df_view.itertuples())

    page = 1
    for i in range(0, len(items), rows_per_page):
        chunk = items[i:i+rows_per_page]
        page_header(c, f"Ligand Summary — page {math.ceil((i+1)/rows_per_page)}")

        table_rows = [["Ligand", "Structure", "best_2A (kcal/mol)", "best_2B (kcal/mol)", "Δ(2A–2B)", "cnn_2A", "cnn_2B"]]
        images = []
        for row in chunk:
            lig = row.Index
            img_path = None
            if lig in lig_mols:
                img = rd_img_from_mol(lig_mols[lig], size=(220,180))
                if img:
                    tmp = Path(tempfile.gettempdir())/f"{lig}.png"
                    img.save(tmp)
                    img_path = str(tmp)
            images.append(img_path)
            table_rows.append([lig, "", f"{row.best_2A:.2f}", f"{row.best_2B:.2f}" if pd.notna(row.best_2B) else "NA",
                               f"{row.delta_2A_2B:.2f}" if pd.notna(row.delta_2A_2B) else "NA",
                               f"{row.cnn_2A:.3f}" if pd.notna(row.cnn_2A) else "NA",
                               f"{row.cnn_2B:.3f}" if pd.notna(row.cnn_2B) else "NA"])

        # build table; we’ll place images after drawing
        colw = [5*cm, 8*cm, 4*cm, 4*cm, 3.5*cm, 3.5*cm, 3.5*cm]
        t = Table(table_rows, colWidths=colw, rowHeights=[0.9*cm] + [3.2*cm]*len(chunk))
        t.setStyle(TableStyle([
            ("FONT", (0,0), (-1,-1), "Helvetica", 9),
            ("FONT", (0,0), (-1,0), "Helvetica-Bold", 10),
            ("BACKGROUND", (0,0), (-1,0), colors.whitesmoke),
            ("GRID", (0,0), (-1,-1), 0.25, colors.grey),
            ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ]))
        w,h = t.wrapOn(c, 28*cm, 14*cm)
        x0, y0 = 2*cm, 17*cm - h
        t.drawOn(c, x0, y0)

        # overlay images in the “Structure” column
        for j, img_path in enumerate(images, start=1):
            if img_path and Path(img_path).exists():
                x = x0 + colw[0] + 0.2*cm
                y = y0 + h - j*3.2*cm + 0.3*cm
                c.drawImage(img_path, x, y, width=colw[1]-0.4*cm, height=2.6*cm, preserveAspectRatio=True)

        c.showPage()

    c.save()

def main():
    import argparse, tempfile
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default=DEFAULT_TAG)
    ap.add_argument("--size", type=int, default=24)
    ap.add_argument("--exh", type=int, default=128)
    ap.add_argument("--modes", type=int, default=30)
    args = ap.parse_args()

    raw, summ = load_tables(args.tag)
    grids = json.loads((REPO/"data/grids/grids.json").read_text())
    lig_mols = load_ligands(args.tag)

    # order ligands by selectivity then potency
    view = summ.copy()
    view["rank_key"] = list(zip(view["delta_2A_2B"], -view["best_2A"]))
    view = view.sort_values(["delta_2A_2B","best_2A"])
    # PDF path
    out_pdf = REPO / f"data/screens/{args.tag}/run_report_{args.tag}.pdf"

    # params page
    add_params_page(out_pdf, args.tag, grids, vars(args))
    # append ligand table pages
    add_table_pages_append(out_pdf, view, lig_mols, vars(args))

    print(f"[done] wrote {out_pdf}")

if __name__ == "__main__":
    main()
