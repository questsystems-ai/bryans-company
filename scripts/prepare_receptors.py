from pathlib import Path
import subprocess, shutil, sys

receptor_dir = Path("data/receptors")
pdbs = sorted(receptor_dir.glob("*.pdb"))
if not pdbs:
    print("[warn] No .pdb receptor files found in data/receptors/")
    sys.exit(0)

# Possible ADFR locations on Windows
candidate_bins = [
    r"C:\Program Files (x86)\ADFRsuite-1.0\bin",
    r"C:\Program Files\ADFRsuite-1.0\bin",
]
adt_cmd = None
for base in candidate_bins:
    bat = Path(base) / "prepare_receptor.bat"
    py4 = Path(base) / "prepare_receptor4.py"
    if bat.exists():
        adt_cmd = str(bat)
        break
    if py4.exists():
        adt_cmd = str(py4)
        break

# If not found in typical dirs, try PATH
if not adt_cmd:
    adt_cmd = shutil.which("prepare_receptor") or shutil.which("prepare_receptor.bat") \
              or shutil.which("prepare_receptor4.py")

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
        subprocess.run(cmd, check=True, shell=True)  # shell=True helps with .bat on Windows
    elif obabel_cmd:
        cmd = [obabel_cmd, str(pdb), "-O", str(pdbqt), "--addhydrogens"]
        print(f"[run-obabel] {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
    else:
        print("[error] Neither ADFR prepare_receptor nor obabel found.")
        sys.exit(1)

print("[done] Receptor prep complete")
