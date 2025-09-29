#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# Ensure inputs exist
test -f data/grids/grids.json || { echo "Missing data/grids/grids.json"; exit 1; }
test -f data/ligands/lib.sdf  || { echo "Missing data/ligands/lib.sdf";  exit 1; }
for r in data/receptors/5HT2A_active_7WC5_clean.pdbqt \
         data/receptors/5HT2A_active_6WGT_clean.pdbqt \
         data/receptors/5HT2B_4IB4_clean.pdbqt; do
  test -f "$r" || { echo "Missing receptor $r"; exit 1; }
done

# jq helpers
cx() { jq -r ".\"$1\".center[0]" data/grids/grids.json; }
cy() { jq -r ".\"$1\".center[1]" data/grids/grids.json; }
cz() { jq -r ".\"$1\".center[2]" data/grids/grids.json; }
sx() { jq -r ".\"$1\".size[0]"   data/grids/grids.json; }
sy() { jq -r ".\"$1\".size[1]"   data/grids/grids.json; }
sz() { jq -r ".\"$1\".size[2]"   data/grids/grids.json; }

run_one () {
  local tag="$1" pdbqt="$2"
  local out="data/results/dock_${tag}"
  mkdir -p "$out"
  gnina \
    --receptor "$pdbqt" \
    --ligand data/ligands/bench.sdf \
    --center_x "$(cx "$tag")" --center_y "$(cy "$tag")" --center_z "$(cz "$tag")" \
    --size_x   "$(sx "$tag")" --size_y   "$(sy "$tag")" --size_z   "$(sz "$tag")" \
    --exhaustiveness 16 --num_modes 20 \
    --cnn_scoring rescore --seed 42 \
    --out "${out}/poses.sdf" \
    --log "${out}/gnina.log"
}

run_one 5HT2A_active_7WC5 data/receptors/5HT2A_active_7WC5_clean.pdbqt
run_one 5HT2A_active_6WGT data/receptors/5HT2A_active_6WGT_clean.pdbqt
run_one 5HT2B_4IB4        data/receptors/5HT2B_4IB4_clean.pdbqt

echo "[done] Docking complete → data/results/"
