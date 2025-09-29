#!/usr/bin/env python3
import os, json, requests, yaml
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load(open(BASE/"config"/"targets.yaml"))
outdir = BASE/"data"/"receptors"
outdir.mkdir(parents=True, exist_ok=True)

def fetch_pdb(pdb_id: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text

meta = []
for rec in cfg["receptors"]:
    pdb_id = rec["pdb_id"]
    name = rec["name"]
    print(f"[fetch] {pdb_id} -> {name}.pdb")
    pdb_txt = fetch_pdb(pdb_id)
    pdb_path = outdir/f"{name}.pdb"
    pdb_path.write_text(pdb_txt)
    entry = dict(rec)
    entry["pdb_path"] = str(pdb_path.relative_to(BASE))
    meta.append(entry)

(open(outdir/"metadata.json","w")).write(json.dumps(meta, indent=2))
print("[done] metadata -> data/receptors/metadata.json")
