#!/usr/bin/env python3
"""
Align every ligand in an input SDF to an LSD co-crystal reference pose,
so docking can run with --local_only around a sane starting orientation.

Usage:
  python scripts/align_to_lsd.py \
      --ref data/custom_ligands/LSD_from_5TVN.pdb \
      --inp data/ligands/lib_controls.sdf \
      --out data/ligands/lib_controls_seeded.sdf \
      --outdir data/ligands/seeded  \
      --skip-names LSD,7LD

Notes:
- Tries MCS/core alignment first, then falls back to shape-based O3A.
- Preserves 3D coords; embeds + MMFF minimization if a molecule lacks conformers.
- Sets props: _Name (kept/normalized), SMILES (computed if missing), ALIGN_METHOD, ALIGN_STATUS.
"""

from __future__ import annotations
import argparse, sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS

def load_ref(path: Path) -> Chem.Mol:
    ext = path.suffix.lower()
    if ext in [".pdb", ".mol", ".mol2"]:
        ref = Chem.MolFromPDBFile(str(path), removeHs=False) if ext==".pdb" else Chem.MolFromMolFile(str(path), removeHs=False)
    else:
        # SDF: take first mol
        sup = Chem.SDMolSupplier(str(path), removeHs=False)
        ref = next((m for m in sup if m), None)
    if not ref:
        sys.exit(f"[ERR] Could not read reference from {path}")
    # Ensure H + coords
    ref = Chem.AddHs(ref, addCoords=True)
    if ref.GetNumConformers()==0:
        AllChem.EmbedMolecule(ref, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(ref)
    return ref

def get_name(m: Chem.Mol, fallback: str) -> str:
    for key in ("_Name","NAME","Name","name"):
        if m.HasProp(key):
            v = m.GetProp(key).strip()
            if v: return v
    return fallback

def ensure_smiles(m: Chem.Mol):
    if not m.HasProp("SMILES"):
        try:
            m.SetProp("SMILES", Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True))
        except Exception:
            pass

def align_to_ref(m: Chem.Mol, ref: Chem.Mol) -> tuple[Chem.Mol,str,bool]:
    """Return aligned copy, method label, and success flag."""
    mm = Chem.Mol(m)  # copy
    # Ensure H + embedding if needed
    mm = Chem.AddHs(mm, addCoords=True)
    if mm.GetNumConformers()==0:
        AllChem.EmbedMolecule(mm, AllChem.ETKDG())
        try:
            AllChem.MMFFOptimizeMolecule(mm)
        except Exception:
            pass
    # Try MCS core alignment
    try:
        res = rdFMCS.FindMCS([ref, mm], timeout=5, ringMatchesRingOnly=True, completeRingsOnly=True)
        if res.smartsString:
            patt = Chem.MolFromSmarts(res.smartsString)
            rt = ref.GetSubstructMatch(patt); mt = mm.GetSubstructMatch(patt)
            if rt and mt:
                AllChem.AlignMol(mm, ref, atomMap=list(zip(mt, rt)))
                return mm, "MCS", True
    except Exception:
        pass
    # Fallback: shape alignment
    try:
        o3a = AllChem.GetO3A(mm, ref)
        o3a.Align()
        return mm, "O3A", True
    except Exception:
        return Chem.Mol(m), "NONE", False

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--ref", required=True, type=Path, help="Reference LSD pose (PDB/SDF/MOL)")
    ap.add_argument("--inp", required=True, type=Path, help="Input SDF with ligands to align")
    ap.add_argument("--out", required=True, type=Path, help="Output seeded SDF")
    ap.add_argument("--outdir", type=Path, default=None, help="Also write one aligned SDF per ligand here")
    ap.add_argument("--skip-names", type=str, default="LSD,7LD", help="Comma-separated names to skip (already reference)")
    args = ap.parse_args()

    ref = load_ref(args.ref)
    skip = {x.strip().lower() for x in args.skip_names.split(",")} if args.skip_names else set()

    suppl = Chem.SDMolSupplier(str(args.inp), removeHs=False)
    W = Chem.SDWriter(str(args.out))
    per_dir = args.outdir
    if per_dir: per_dir.mkdir(parents=True, exist_ok=True)

    count=good=0
    for i, m in enumerate(suppl, start=1):
        if m is None: 
            continue
        count += 1
        name = get_name(m, f"mol_{i}")
        m.SetProp("_Name", name)
        ensure_smiles(m)
        # Skip if this is the reference itself
        if name.lower() in skip:
            m.SetProp("ALIGN_METHOD","SKIP")
            m.SetProp("ALIGN_STATUS","SKIPPED")
            W.write(m)
            if per_dir: Chem.SDWriter(str(per_dir/f"{name}.sdf")).write(m)
            continue
        aligned, method, ok = align_to_ref(m, ref)
        aligned.SetProp("_Name", name)
        ensure_smiles(aligned)
        aligned.SetProp("ALIGN_METHOD", method)
        aligned.SetProp("ALIGN_STATUS", "OK" if ok else "FAIL")
        W.write(aligned)
        if per_dir:
            Chem.SDWriter(str(per_dir/f"{name}.sdf")).write(aligned)
        if ok: good += 1

    W.close()
    print(f"[done] wrote {args.out}  ({good}/{count} aligned OK)")
    if per_dir:
        print(f"[note] individual aligned SDFs in: {per_dir}")

if __name__ == "__main__":
    main()
