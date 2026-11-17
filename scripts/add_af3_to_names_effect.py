"""Add AlphaFold3 to data/names_effect_per_structure.csv, in the same format as the other three models.

    python scripts/add_af3_to_names_effect.py        (safe to re-run: it replaces any existing AlphaFold3 rows)

For each structure: did the best of the 5 samples at seed 42 land within 2 A, under each input. Two steps, like the others:
    chemistry only             smiles -> gen   (gen = the CCD molecule with generic atom names, AlphaFold3's ccd_genericnames)
    names on top of chemistry  gen -> ccd      (generic names -> deposited names)
Source: data/pxr_cyp3a4_pose_rmsds.csv.gz (AlphaFold3's poses for all three inputs).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import common as C

SEED = 42
FILE = C.DATA / "names_effect_per_structure.csv"
ARM = {"smiles": "smiles", "gen": "ccd_genericnames", "ccd": "ccd"}      # name used in the file -> input in the poses file
STEPS = [("chemistry only", "smiles", "gen"), ("names on top of chemistry", "gen", "ccd")]

poses = pd.read_csv(C.DATA / "pxr_cyp3a4_pose_rmsds.csv.gz")
rows = []
for tgt_label, tgt in (("PXR", "pxr"), ("CYP3A4", "cyp3a4")):
    a = poses[(poses.model == "af3") & (poses.target == tgt) & (poses.seed == SEED) & (poses.status == "ok")]
    best = a.groupby(["structure", "input"]).ligand_rmsd.min().unstack("input")
    hit = (best <= 2)
    for step, lo, hi in STEPS:
        pair = best[[ARM[lo], ARM[hi]]].dropna()           # structures that ran in both arms of this step
        for s in sorted(pair.index):
            x, y = int(hit.loc[s, ARM[lo]]), int(hit.loc[s, ARM[hi]])
            rows.append(dict(target=tgt_label, model="AlphaFold3", step=step, structure=s, lower_arm=lo, upper_arm=hi,
                             lower_hit=x, upper_hit=y, direction="improved" if y > x else ("worsened" if y < x else "unchanged")))

old = pd.read_csv(FILE)
old = old[old.model != "AlphaFold3"]
out = pd.concat([old, pd.DataFrame(rows)], ignore_index=True)
out.to_csv(FILE, index=False, lineterminator="\n")
print(out.groupby(["model", "target", "step"]).size().to_string())
