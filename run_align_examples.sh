#!/usr/bin/env bash
set -euo pipefail

REPO=~/5ht-platform-starter
REF="$REPO/data/custom_ligands/LSD_from_5TVN.pdb"   # or LSD.sdf if you saved it
IN1="$REPO/data/ligands/lib_controls.sdf"
OUT1="$REPO/data/ligands/lib_controls_seeded.sdf"

IN2="$REPO/data/ligands/controls_trio_pH74.sdf"
OUT2="$REPO/data/ligands/controls_trio_pH74_seeded.sdf"

mkdir -p "$REPO/data/ligands/seeded"

python "$REPO/scripts/align_to_lsd.py" \
  --ref "$REF" \
  --inp "$IN1" \
  --out "$OUT1" \
  --outdir "$REPO/data/ligands/seeded" \
  --skip-names LSD,7LD

python "$REPO/scripts/align_to_lsd.py" \
  --ref "$REF" \
  --inp "$IN2" \
  --out "$OUT2" \
  --outdir "$REPO/data/ligands/seeded" \
  --skip-names LSD,7LD

echo "Aligned SDFs:"
ls -lh "$REPO/data/ligands/" | grep seeded || true
