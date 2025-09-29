#!/usr/bin/env bash
set -euo pipefail

BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
RECP="$BASE/data/receptors"
GRID="$BASE/data/grids/grids.json"
LIGS="$BASE/data/ligands/lib.sdf"
OUTD="$BASE/data/results"

if [ ! -f "$GRID" ]; then
  echo "[err] $GRID not found. Run scripts/make_grids_from_cocrystal.py first."
  exit 1
fi
if [ ! -f "$LIGS" ]; then
  echo "[err] $LIGS not found. Run scripts/prepare_ligands.py first."
  exit 1
fi

EXH=16
MODES=20

docker run --rm -u $(id -u):$(id -g) \
  -v "$BASE":"$BASE" \
  gnina/gnina:latest bash -lc '
python - <<PY
import json, os, subprocess, pathlib
BASE = pathlib.Path("'"$BASE"'")
grids = json.load(open(BASE/"data"/"grids"/"grids.json"))
for name, g in grids.items():
    rec_pdbqt = BASE/"data"/"receptors"/(pathlib.Path(g["pdb"]).stem + ".pdbqt")
    if not rec_pdbqt.exists():
        print(f"[warn] {rec_pdbqt.name} missing; skip {name}")
        continue
    outdir = BASE/"data"/"results"/f"dock_{name}"
    outdir.mkdir(parents=True, exist_ok=True)
    cx, cy, cz = g["center"]
    sx, sy, sz = g["size"]
    cmd = [
        "gnina",
        "--receptor", str(rec_pdbqt),
        "--ligand", str(BASE/"data"/"ligands"/"lib.sdf"),
        "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
        "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
        "--exhaustiveness", str('"$EXH"'),
        "--num_modes", str('"$MODES"'),
        "--cnn_scoring", "rescore",
        "--seed", "42",
        "--out", str(outdir/"poses.sdf"),
        "--log", str(outdir/"gnina.log")
    ]
    print("[dock]", name)
    subprocess.run(cmd, check=True)
PY
'
echo "[done] docking complete -> $OUTD"
