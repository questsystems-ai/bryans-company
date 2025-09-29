# scripts/prepare_receptors.py
from pathlib import Path
import subprocess, shutil, sys

receptor_dir = Path("data/receptors")

# Prefer *_clean.pdb if present, else fallback to *.pdb
pdbs = sorted(receptor_dir.glob("*_clean.pdb")) or sorted(receptor_dir.glob("*.pdb"))
if not pdbs:
    print("[warn] No receptor .pdb files found in data/receptors/")
    sys.exit(0)

# Locate ADFR (Windows uses .bat) or fallback to obabel
adt_cmd = (
    shutil.which("prepare_receptor")
    or shutil.which("prepare_receptor.bat")
    or shutil.which("prepare_receptor4.py")
)
obabel_cmd = shutil.which("obabel")

print("[info] ADFR prep:", adt_cmd or "not found")
print("[info] obabel   :", obabel_cmd or "not found")

for pdb in pdbs:
    pdbqt = pdb.with_suffix(".pdbqt")
    if pdbqt.exists():
        print(f"[skip] {pdbqt} already exists")
        continue

    if adt_cmd:
        cmd = [adt_cmd, "-r", str(pdb), "-o", str(pdbqt), "-A", "hydrogens"]
        print(f"[run-adt] {' '.join(cmd)}")
        subprocess.run(cmd, check=True, shell=True)  # shell=True helps with .bat
    elif obabel_cmd:
        cmd = [obabel_cmd, str(pdb), "-O", str(pdbqt), "--addhydrogens"]
        print(f"[run-obabel] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    else:
        print("[error] Neither prepare_receptor nor obabel found.")
        sys.exit(1)

print("[done] Receptor prep complete")
