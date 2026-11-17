"""Write data/structure_ids.csv: every PDB entry behind every figure, with the figures it appears in.

    python scripts/make_structure_ids.py        (run from anywhere; needs only the CSVs in data/)

Columns: set, pdb_id, deposit_date, release_date, ligand_ccd (random sets only), figures.
  pxr / cyp3a4   the common sets (49 + 85): every model and input produced poses on all five seeds (Figs 1, 4)
                 "3" is added for those also in Figure 3 (shared by all four models in all three arms, seed 42)
  random_before / random_after   the 50 + 50 random complexes in Figure 5
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import common as C

pdb = lambda s: s.split("_")[1].upper()
B = C.load_poses()
dates = pd.read_csv(C.DATA / "matrix_dates.csv").set_index("pdb")

# Figure 3: structures every model ran in all three arms at seed 42 (same rule as the Figure 3 notebook)
ps = pd.read_csv(C.DATA / "names_effect_per_structure.csv")
poses = pd.read_csv(C.DATA / "matrix_poses.csv.gz")
fig3 = {}
for tgt, key in (("PXR", "pxr"), ("CYP3A4", "cyp3a4")):
    sets = []
    for m in ("Boltz-2", "OF3 p2", "OBv0"):
        s = ps[(ps.model == m) & (ps.target == tgt)]
        sets.append(set(s[s.step == "chemistry only"].structure) & set(s[s.step == "names on top of chemistry"].structure))
    a = poses[(poses.model == "af3") & (poses.status == "ok") & (poses.target == key) & (poses.seed == 42)]
    best = a.groupby(["input", "structure"]).ligand_rmsd.min().unstack("input")[["smiles", "ccd_genericnames", "ccd"]].dropna()
    sets.append(set(best.index))
    fig3[key] = set.intersection(*sets)
assert (len(fig3["pxr"]), len(fig3["cyp3a4"])) == (48, 85), {k: len(v) for k, v in fig3.items()}

rows = []
for t in ("pxr", "cyp3a4"):
    for s in sorted(B[B.target == t].structure.unique()):
        p = pdb(s)
        rows.append(dict(set=t, pdb_id=p, deposit_date=dates.loc[p, "deposit"], release_date=dates.loc[p, "release"],
                         ligand_ccd="", figures="1,3,4" if s in fig3[t] else "1,4"))

info = pd.read_csv(C.DATA / "random_complexes_info.csv").set_index("pdb_id")
for side in ("before", "after"):
    ids = sorted(pd.read_csv(C.DATA / f"random_complexes_{side}_cutoff.csv").pdb_id.unique())
    assert len(ids) == 50
    for p in ids:
        rows.append(dict(set=f"random_{side}", pdb_id=p, deposit_date=info.loc[p, "deposit_date"], release_date=info.loc[p, "release_date"],
                         ligand_ccd=info.loc[p, "ligand_ccd"], figures="5"))

out = pd.DataFrame(rows)
out.to_csv(C.DATA / "structure_ids.csv", index=False)
print(out.groupby("set").size().to_dict())
