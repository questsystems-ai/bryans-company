#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO"

CPU_THREADS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpu) CPU_THREADS="$2"; shift 2;;
    *) echo "Unknown arg: $1"; exit 2;;
  esac
done

[[ -f data/grids/grids.json ]] || { echo "Missing data/grids/grids.json"; exit 1; }
[[ -f data/ligands/bench.sdf ]] || { echo "Missing data/ligands/bench.sdf"; exit 1; }
for r in data/receptors/5HT2A_active_7WC5_clean.pdbqt \
         data/receptors/5HT2A_active_6WGT_clean.pdbqt \
         data/receptors/5HT2B_4IB4_clean.pdbqt; do
  [[ -f "$r" ]] || { echo "Missing receptor $r"; exit 1; }
done

if ! command -v jq >/dev/null 2>&1; then
  echo "Installing jq..."
  sudo apt-get update -y && sudo apt-get install -y jq
fi

cx(){ jq -r ".\"$1\".center[0]" data/grids/grids.json; }
cy(){ jq -r ".\"$1\".center[1]" data/grids/grids.json; }
cz(){ jq -r ".\"$1\".center[2]" data/grids/grids.json; }
sx(){ jq -r ".\"$1\".size[0]"   data/grids/grids.json; }
sy(){ jq -r ".\"$1\".size[1]"   data/grids/grids.json; }
sz(){ jq -r ".\"$1\".size[2]"   data/grids/grids.json; }

DOCKER_GPU_FLAGS="--gpus all"
EXTRA_GNINA_FLAGS=""
if [[ -n "${CPU_THREADS}" ]]; then
  DOCKER_GPU_FLAGS=""
  EXTRA_GNINA_FLAGS="--cpu ${CPU_THREADS}"
fi

run_one() {
  local tag="$1" pdbqt="$2"
  local out="data/results/dock_${tag}"
  mkdir -p "$out"
  local CX="$(cx "$tag")" CY="$(cy "$tag")" CZ="$(cz "$tag")"
  local SX="$(sx "$tag")" SY="$(sy "$tag")" SZ="$(sz "$tag")"
  echo ">>> Docking ${tag}  [center=(${CX},${CY},${CZ}) size=(${SX},${SY},${SZ})]"
  docker run --rm -it ${DOCKER_GPU_FLAGS} \
    -v "$REPO":/repo \
    gnina/gnina:latest \
    gnina \
      --receptor /repo/"$pdbqt" \
      --ligand   /repo/data/ligands/bench.sdf \
      --center_x "${CX}" --center_y "${CY}" --center_z "${CZ}" \
      --size_x   "${SX}" --size_y   "${SY}" --size_z   "${SZ}" \
      --exhaustiveness 16 --num_modes 20 \
      --cnn_scoring rescore --seed 42 \
      ${EXTRA_GNINA_FLAGS} \
      --out /repo/"${out}/poses.sdf" \
      --log /repo/"${out}/gnina.log"
}

run_one 5HT2A_active_7WC5 data/receptors/5HT2A_active_7WC5_clean.pdbqt
run_one 5HT2A_active_6WGT data/receptors/5HT2A_active_6WGT_clean.pdbqt
run_one 5HT2B_4IB4        data/receptors/5HT2B_4IB4_clean.pdbqt

SUMMARY="data/results_summary.csv"
echo "target,pose_rank,affinity,cnnscore,minimizedAffinity,file" > "$SUMMARY"
summarize() {
  local sdf="$1" tag="$2"
  awk -v T="$tag" -v F="$sdf" '
    /^REMARK Vina result:/ { match($0,/result:\s*([-0-9.]+)/,m); aff=m[1]; next }
    /^>  <CNNscore>/ { getline; cnn=$0; next }
    /^>  <minimizedAffinity>/ { getline; minA=$0; next }
    /^$$$$/ { pose++; printf "%s,%d,%s,%s,%s,%s\n", T, pose, aff, cnn, minA, F; aff=cnn=minA="" }
  ' "$sdf"
}
for d in data/results/dock_*; do
  [[ -f "$d/poses.sdf" ]] && summarize "$d/poses.sdf" "${d##*/dock_}" >> "$SUMMARY"
done

echo "[done] Results in data/results/ ; summary: ${SUMMARY}"
