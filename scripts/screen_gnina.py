#!/usr/bin/env python3
import argparse, json, re, subprocess, sys, tempfile, shutil
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

RX_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d\.\d+)\s+(-?\d+\.\d+)")

TARGETS = ["5HT2A_active_7WC5", "5HT2A_active_6WGT", "5HT2B_4IB4"]

def run(cmd, **kw):
    p = subprocess.run(cmd, shell=isinstance(cmd,str), capture_output=True, text=True, **kw)
    if p.returncode != 0:
        print(p.stdout, p.stderr, file=sys.stderr)
        raise RuntimeError(f"cmd failed: {cmd}")
    return p

def prep_protonated_sdf(seed_smi: Path, out_sdf: Path):
    # Use Dimorphite 1.x CLI to protonate, then RDKit to 3D
    prot_smi = out_sdf.with_suffix(".prot.smi")
    run(["dimorphite_dl", "--min_ph", "6.8", "--max_ph", "8.0",
         "--smiles_file", str(seed_smi), "--output_file", str(prot_smi), "--silent"])
    w = Chem.SDWriter(str(out_sdf)); n=0
    for line in prot_smi.read_text().splitlines():
        s=line.strip()
        if not s: continue
        parts=s.split()
        smi = parts[0]
        name = parts[1] if len(parts)>1 else f"lig_{n+1}"
        m = Chem.MolFromSmiles(smi)
        if not m: continue
        m = Chem.AddHs(m)
        try:
            AllChem.EmbedMolecule(m, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass
        m.SetProp("_Name", name)
        w.write(m); n+=1
    w.close()
    if n == 0:
        raise RuntimeError("no molecules written to SDF")
    return n

def split_sdf(multi: Path, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    suppl = Chem.SDMolSupplier(str(multi), removeHs=False)
    paths=[]
    for i,m in enumerate(suppl):
        if m is None: continue
        name = m.GetProp("_Name") if m.HasProp("_Name") else f"lig_{i+1}"
        out = out_dir / f"{name}.sdf"
        w = Chem.SDWriter(str(out))
        w.write(m); w.close()
        paths.append(out)
    if not paths:
        raise RuntimeError("split produced no SDFs")
    return paths

def parse_gnina_log(log_path: Path):
    best = None
    # returns list of rows per mode
    out=[]
    with log_path.open(errors="ignore") as fh:
        in_table=False
        for line in fh:
            if line.strip().startswith("mode |") and "affinity" in line:
                in_table=True; continue
            if in_table:
                m = RX_ROW.match(line)
                if m:
                    mode      = int(m.group(1))
                    affinity  = float(m.group(2))
                    intramol  = float(m.group(3))
                    cnn_pose  = float(m.group(4))
                    cnn_aff   = float(m.group(5))
                    out.append(dict(mode=mode, affinity=affinity, intramol=intramol,
                                    cnnpose=cnn_pose, cnnaff=cnn_aff))
                elif line.strip()=="" and out:
                    in_table=False
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sdf", default="data/ligands/lib.sdf", help="input ligands SDF")
    ap.add_argument("--smi", default="data/ligands/seed.smi", help="if SDF missing, protonate this .smi to SDF")
    ap.add_argument("--img", default="gnina/gnina:latest", help="docker image")
    ap.add_argument("--exh", type=int, default=128)
    ap.add_argument("--modes", type=int, default=30)
    ap.add_argument("--tag", default="screen1", help="output tag under data/screens/<tag>")
    ap.add_argument("--size", type=int, default=24, help="box size edge (Å)")
    args = ap.parse_args()

    repo = Path.cwd()
    sdf = repo/args.sdf
    if not sdf.exists():
        print(f"[info] {sdf} not found; protonating {args.smi} -> {sdf}")
        prep_protonated_sdf(repo/args.smi, sdf)

    # split to one-ligand files
    lig_dir = repo / f"data/screens/{args.tag}/ligands"
    lig_paths = split_sdf(sdf, lig_dir)

    # load grids
    grids = json.loads((repo/"data/grids/grids.json").read_text())

    rows=[]
    print(f"[run] ligands: {len(lig_paths)} | exh={args.exh} modes={args.modes} size={args.size}")
    for lig in tqdm(lig_paths, ncols=80):
        lig_name = lig.stem
        for tgt in TARGETS:
            cx,cy,cz = grids[tgt]["center"]
            outdir = repo / f"data/screens/{args.tag}/dock_{tgt}/{lig_name}"
            outdir.mkdir(parents=True, exist_ok=True)
            cmd = [
                "docker","run","--rm","--gpus","all",
                "-v", f"{repo}:{repo}", "-w", f"{repo}",
                args.img, "gnina",
                "--receptor", f"data/receptors/{tgt}_clean.pdbqt",
                "--ligand",   str(lig.relative_to(repo)),
                "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                "--size_x", str(args.size), "--size_y", str(args.size), "--size_z", str(args.size),
                "--exhaustiveness", str(args.exh),
                "--num_modes", str(args.modes),
                "--cnn_scoring", "rescore",
                "--out", str((outdir/"poses.sdf").relative_to(repo)),
                "--log", str((outdir/"gnina.log").relative_to(repo)),
            ]
            run(cmd)

            # parse log
            modes = parse_gnina_log(outdir/"gnina.log")
            if modes:
                best = sorted(modes, key=lambda r: r["affinity"])[0]
                rows.append(dict(
                    ligand=lig_name, target=tgt,
                    affinity=best["affinity"], cnnpose=best["cnnpose"], cnnaff=best["cnnaff"],
                    log=str((outdir/"gnina.log").relative_to(repo))
                ))
            else:
                rows.append(dict(ligand=lig_name, target=tgt, affinity=float("nan"),
                                 cnnpose=float("nan"), cnnaff=float("nan"),
                                 log=str((outdir/"gnina.log").relative_to(repo))))

    df = pd.DataFrame(rows)
    # compute best 2A per ligand (min across 7WC5/6WGT)
    pivot = df.pivot_table(index="ligand", columns="target", values="affinity", aggfunc="min")
    df2 = pivot.copy()
    df2["best_2A"] = df2[["5HT2A_active_7WC5","5HT2A_active_6WGT"]].min(axis=1)
    df2["best_2B"] = df2["5HT2B_4IB4"]
    df2["delta_2A_2B"] = df2["best_2A"] - df2["best_2B"]  # more negative = more selective for 2A
    df2 = df2.sort_values("delta_2A_2B")  # most selective on top

    out_dir = repo / f"data/screens/{args.tag}"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir/"raw_best_by_target.csv", index=False)
    df2.to_csv(out_dir/"summary_per_ligand.csv")  # index=ligand
    print(f"[done] wrote:\n  {out_dir/'raw_best_by_target.csv'}\n  {out_dir/'summary_per_ligand.csv'}")

if __name__ == "__main__":
    main()
