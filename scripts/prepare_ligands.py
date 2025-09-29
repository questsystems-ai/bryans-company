# scripts/prepare_ligands.py
import argparse
import subprocess
import sys
import shutil
from pathlib import Path
from typing import List

from rdkit import Chem
from rdkit.Chem import AllChem


# --- Dimorphite-DL compatibility (v1.x and v2.x) ---
import shutil, subprocess, sys
from typing import List

def _import_dimorphite():
    """
    Return (api_tag, DimorphiteDL_cls or None).
    Try nested (v1-style) first, then top-level (some v2 builds).
    """
    try:
        from dimorphite_dl.dimorphite_dl import DimorphiteDL  # v1 path (often also present)
        return ("v1", DimorphiteDL)
    except Exception:
        pass
    try:
        from dimorphite_dl import DimorphiteDL  # v2 top-level (may or may not export)
        return ("v2", DimorphiteDL)
    except Exception:
        return (None, None)

def _protonate_via_import(smi: str, min_ph: float, max_ph: float, max_variants: int) -> List[str]:
    api, DimorphiteDL = _import_dimorphite()
    if DimorphiteDL is None:
        return []
    engine = DimorphiteDL(min_ph=min_ph, max_ph=max_ph, pka_precision=0.5)
    try:
        # v2 signature supports kwargs + max_variants
        variants = engine.protonate(smiles=smi, min_ph=min_ph, max_ph=max_ph, max_variants=max_variants)
    except TypeError:
        # v1 signature is positional only
        variants = engine.protonate(smi)
    return list(dict.fromkeys(variants))

def _protonate_via_cli(smi: str, min_ph: float, max_ph: float, max_variants: int) -> List[str]:
    """
    Use the console script 'dimorphite_dl' if imports aren't available.
    Tries v2 flags first, then v1 flags.
    """
    exe = shutil.which("dimorphite_dl")
    if not exe:
        return []
    # --- Try v2.x style ---
    try:
        cmd2 = [
            exe,
            "--ph_min", str(min_ph),
            "--ph_max", str(max_ph),
            "--precision", "0.5",
            "--max_variants", str(max_variants),
            smi,  # positional SMI
        ]
        out = subprocess.check_output(cmd2, text=True).strip()
        v = [line.split()[0] for line in out.splitlines() if line.strip()]
        if v:
            return list(dict.fromkeys(v))
    except subprocess.CalledProcessError:
        pass
    # --- Fallback to v1.x style ---
    try:
        cmd1 = [
            exe,
            "--min_ph", str(min_ph),
            "--max_ph", str(max_ph),
            "--pka_precision", "0.5",
            "--smiles", smi,
            "--output_format", "smi",
        ]
        out = subprocess.check_output(cmd1, text=True).strip()
        v = [line.split()[0] for line in out.splitlines() if line.strip()]
        if v:
            return list(dict.fromkeys(v))
    except subprocess.CalledProcessError:
        pass
    return []

def protonate_smiles(smi: str, min_ph: float, max_ph: float, max_variants: int) -> List[str]:
    # Prefer Python API; if not available, use CLI (v2 then v1)
    v = _protonate_via_import(smi, min_ph, max_ph, max_variants)
    if v:
        return v
    v = _protonate_via_cli(smi, min_ph, max_ph, max_variants)
    if v:
        return v
    raise RuntimeError("Dimorphite-DL not available (import and CLI both failed).")

def write_sdf(mols: List[Chem.Mol], out_path: Path):
    w = Chem.SDWriter(str(out_path))
    for m in mols:
        w.write(m)
    w.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smi", required=True, help="Input SMILES file (one per line; optional name after SMILES)")
    ap.add_argument("--out", required=True, help="Output SDF path")
    ap.add_argument("--num_confs", type=int, default=20)
    ap.add_argument("--min_ph", type=float, default=6.8)
    ap.add_argument("--max_ph", type=float, default=8.0)
    ap.add_argument("--max_variants", type=int, default=128)
    args = ap.parse_args()

    in_path = Path(args.smi)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    out_mols: List[Chem.Mol] = []

    with in_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            base = parts[0]
            name = parts[1] if len(parts) > 1 else base

            try:
                variants = protonate_smiles(base, args.min_ph, args.max_ph, args.max_variants)
            except Exception as e:
                print(f"[warn] protonation failed for {base}: {e}", file=sys.stderr)
                continue

            for i, smi in enumerate(variants):
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    continue
                m = Chem.AddHs(m)
                AllChem.EmbedMultipleConfs(
                    m,
                    numConfs=args.num_confs,
                    useExpTorsionAnglePrefs=True,
                    useBasicKnowledge=True,
                )
                AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=200)
                m.SetProp("_Name", f"{name}__p{i}")
                m.SetProp("input_smiles", base)
                m.SetProp("prot_variant_index", str(i))
                out_mols.append(m)

    if not out_mols:
        sys.exit("[error] No valid ligands were produced. Check SMILES and Dimorphite-DL install.")

    write_sdf(out_mols, out_path)
    print(f"[done] Wrote {len(out_mols)} molecules to {out_path}")


if __name__ == "__main__":
    main()
