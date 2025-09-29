#!/usr/bin/env python3
import json, yaml
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parents[1]
meta = json.load(open(BASE/"data"/"receptors"/"metadata.json"))
cfg = yaml.safe_load(open(BASE/"config"/"targets.yaml"))
grid_size = np.array(cfg["grid"]["size"], dtype=float)

def parse_pdb_for_ligand_coords(path):
    coords = []
    protein_res = {"ALA","VAL","LEU","ILE","MET","PHE","TYR","TRP","PRO","SER","THR","CYS","ASN","GLN","ASP","GLU","HIS","LYS","ARG","GLY","SEC","PYL"}
    for line in open(path):
        if line.startswith(("HETATM","ATOM  ")):
            resn = line[17:20].strip()
            if resn in ("HOH","WAT"): 
                continue
            if resn not in protein_res:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                coords.append((x,y,z))
    if not coords:
        raise RuntimeError("No non-protein HETATM atoms found to define grid center.")
    arr = np.array(coords)
    return arr.mean(axis=0)

grids = {}
for entry in meta:
    pdb_path = BASE/entry["pdb_path"]
    center = parse_pdb_for_ligand_coords(pdb_path)
    grids[entry["name"]] = {
        "center": center.tolist(),
        "size": grid_size.tolist(),
        "pdb": entry["pdb_path"],
        "pdb_id": entry["pdb_id"],
        "state": entry["state"]
    }

gpath = BASE/"data"/"grids"/"grids.json"
gpath.write_text(json.dumps(grids, indent=2))
print(f"[done] grids -> {gpath.relative_to(BASE)}")
