"""Shared helpers for every figure notebook: data loading, the two statistics, and the bar-chart style.

Success = the best of the 5 samples in a seed has ligand RMSD < 2 A. A structure's score is the fraction of
its 5 seeds that succeed; a model's success rate is the mean of those scores over structures.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, wilcoxon

DATA = Path(__file__).resolve().parent.parent / "data"
FIGS = Path(__file__).resolve().parent.parent / "figures"
FIGS.mkdir(exist_ok=True)

# colours, one per ligand input
SMILES, CCD, GENERIC = "#2a78d6", "#eb6834", "#1baf7a"
GREY, INK = "#747d86", "#0b0b0b"

MODELS = ("boltz2", "af3", "of3", "obv0")
LABEL = {"boltz2": "Boltz-2", "af3": "AlphaFold3", "of3": "OF3", "obv0": "OBv0"}

# training cutoffs: OF3 p2 shares AlphaFold3's; OBv0's (2025-06-30) is after nearly every structure, so it is not split
CUTOFF = {"boltz2": pd.Timestamp("2023-06-01"), "af3": pd.Timestamp("2021-09-30"),
          "of3": pd.Timestamp("2021-09-30"), "obv0": pd.Timestamp("2025-06-30")}

# the models x inputs that must all have produced valid poses for a structure to enter the common set
_ARMS = [("af3", "smiles"), ("af3", "ccd"), ("af3", "ccd_genericnames"), ("boltz2", "smiles"), ("boltz2", "ccd"),
         ("obv0", "smiles"), ("obv0", "ccd"), ("of3", "smiles"), ("of3", "ccd")]


def load_poses():
    """One row per (model, target, input, seed, structure): the lowest ligand RMSD over the 5 samples, on the common set
    of structures (49 PXR, 85 CYP3A4) that every model and input produced poses for on all 5 seeds. Columns `ok`
    (RMSD < 2 A), `release`."""
    d = pd.read_csv(DATA / "pxr_cyp3a4_pose_rmsds.csv.gz")
    d = d[d.target.isin(["pxr", "cyp3a4"]) & d.input.isin(["smiles", "ccd", "ccd_genericnames"])]
    b = d.groupby(["model", "target", "input", "seed", "structure"]).ligand_rmsd.min().reset_index()
    keep = None
    for m, i in _ARMS:
        g = b[(b.model == m) & (b.input == i)].pivot_table(index=["target", "structure"], columns="seed", values="ligand_rmsd")
        ok = set(g.index[g.notna().all(axis=1)]) if g.shape[1] == 5 else set()
        keep = ok if keep is None else keep & ok
    b = b[[(t, s) in keep for t, s in zip(b.target, b.structure)]].copy()
    b["ok"] = b.ligand_rmsd < 2
    pdb = b.structure.str.split("_").str[1].str.upper()
    dates = pd.read_csv(DATA / "pdb_pxr_cyp3a4_release_dates.csv", parse_dates=["release"]).set_index("pdb")
    b["release"] = pdb.map(dates.release)
    return b


def wide(b, model, target=None):
    """structures x (input, seed) success matrix for one model."""
    x = b[(b.model == model) & b.input.isin(["smiles", "ccd"])]
    if target:
        x = x[x.target == target]
    return x.pivot_table(index=["target", "structure", "release"], columns=["input", "seed"], values="ok")


def pair_stats(x):
    """SMILES vs CCD for one `wide` table: mean success of each, and the paired two-sided Wilcoxon signed-rank p-value
    across structures on their seed-averaged success (structures with identical scores under both inputs are dropped)."""
    s, c = x["smiles"].mean(axis=1), x["ccd"].mean(axis=1)
    p = wilcoxon(s, c).pvalue if (s != c).any() else np.nan
    return dict(n=len(x), smiles=100 * s.mean(), ccd=100 * c.mean(), p=float(p))


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values (same order as the input)."""
    return list(false_discovery_control(list(pvalues), method="bh"))


def side_mask(x, model, side):
    """Rows of a `wide` table released on/before the model's training cutoff (so it could have been seen in training) or after it."""
    c = CUTOFF[model]
    rel = x.index.get_level_values("release")
    return np.asarray(rel <= c) if side == "before" else np.asarray(rel > c)


def grouped_bars(ax, rows, title, subtitle, xlabel="success: best of 5 within 2 Å (%)"):
    """Horizontal SMILES/CCD bar pairs. rows = [(label, smiles %, ccd %), ...]."""
    h, y = 0.36, np.arange(len(rows))
    ax.barh(y - h / 2, [r[1] for r in rows], h, color=SMILES)
    ax.barh(y + h / 2, [r[2] for r in rows], h, color=CCD)
    for i, (_, s, c) in enumerate(rows):
        ax.text(s + 1, i - h / 2, f"{s:.0f}%", va="center", fontsize=10)
        ax.text(c + 1, i + h / 2, f"{c:.0f}%", va="center", fontsize=10)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=10)
    ax.set_xlim(0, 100)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=10 if xlabel.startswith("success") else 9)
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=18)
    ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=10, color=GREY)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
