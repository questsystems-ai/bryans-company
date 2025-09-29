#!/usr/bin/env python3
import os, json
from pathlib import Path
import pandas as pd
from rdkit import Chem

BASE = Path(__file__).resolve().parents[1]
resdir = BASE/"data"/"results"

rows = []
for dockdir in sorted(resdir.glob("dock_*")):
    name = dockdir.name.replace("dock_","")
    sdf = dockdir/"poses.sdf"
    if not sdf.exists():
        continue
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    for i, m in enumerate(suppl):
        if m is None: continue
        lig = m.GetProp("_Name") if m.HasProp("_Name") else f"lig_{i}"
        vina = float(m.GetProp("minimizedAffinity")) if m.HasProp("minimizedAffinity") else float(m.GetProp("affinity")) if m.HasProp("affinity") else None
        cnn = float(m.GetProp("CNNscore")) if m.HasProp("CNNscore") else None
        rows.append({"receptor": name, "ligand": lig, "vina": vina, "cnn": cnn, "rank": i+1})
df = pd.DataFrame(rows)
if df.empty:
    print("[warn] no poses found"); raise SystemExit

df = df.sort_values(["ligand","receptor","cnn","vina"], ascending=[True,True,False,True]).groupby(["ligand","receptor"], as_index=False).first()

df_piv = df.pivot(index="ligand", columns="receptor", values="cnn")
df_piv.to_parquet(BASE/"data"/"results"/"best_cnn_by_receptor.parquet")
df.to_parquet(BASE/"data"/"results"/"best_poses.parquet", index=False)

print("[done] wrote:")
print(" - data/results/best_poses.parquet")
print(" - data/results/best_cnn_by_receptor.parquet")
