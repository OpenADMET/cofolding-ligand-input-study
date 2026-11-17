# cofolding-ligand-input-study

Code, data and figures to accompany the blog post **"If You Give A Cofolding Model a CCD ID"** (Kate Huddleston and Pat Walters).

We gave Boltz-2, AlphaFold3, OpenFold3 (OF3) and OpenFold-OpenBind v0 (OBv0) the same ligand as a SMILES string or as a PDB chemical-component (CCD) code, on PXR, CYP3A4 and a random set of PDB complexes, and compared how often the best predicted pose landed within 2 Å of the experimental ligand.

## Figures

One notebook per figure in the post. Each reads from `data/`, prints the table behind the figure, and writes the PNG to `figures/`.

| figure | notebook | shows |
|---|---|---|
| 1 | `notebooks/fig1_smiles_vs_ccd.ipynb` | SMILES vs CCD success per model and target; BH-adjusted Wilcoxon p-values; per-seed dots |
| 2 | `notebooks/fig2_input_schematic.ipynb` | the same ligand (oxadiazon, CCD S6Z) as SMILES and as a CCD code |
| 3 | `notebooks/fig3_atom_names.ipynb` | what the deposited atom names add (single seed) |
| 4 | `notebooks/fig4_cutoff_split.ipynb` | the SMILES/CCD gap before vs after each model's training cutoff |
| 5 | `notebooks/fig5_random_complexes.ipynb` | the same split on 50 random PDB complexes (single seed) |

The regenerated PNGs in `figures/` are identical to the ones in the post.

## Run it

```bash
mamba env create -f environment.yml
mamba activate cofolding-ligand-input-study
scripts/run_all.sh          # or open any notebook in notebooks/ and run it top to bottom
```

## Layout

```
data/          processed inputs and structure_ids.csv (see data/README.md)
structures/    the two files Figure 2 needs
src/common.py  data loading, the two statistics, bar-chart style (imported by the notebooks)
notebooks/     one notebook per figure
figures/       output PNGs
scripts/       run_all.sh, make_structure_ids.py
```

## Definitions

- **Success:** the lowest-RMSD of the 5 samples in a seed is within 2 Å of the experimental ligand pose.
- **Structure score:** the fraction of its 5 seeds that succeed. **Model success rate:** the mean score over structures.
- **Common set (Figs 1 and 4):** the 49 PXR and 85 CYP3A4 structures for which every model and input produced poses on all five seeds.
- **p-values (Fig 1 only):** paired two-sided Wilcoxon signed-rank test across structures on seed-averaged success, SMILES vs CCD, structures with identical scores dropped, Benjamini-Hochberg adjusted across the eight model-target comparisons.
- **Figs 3 and 5** use one seed (five samples), so their baselines differ from Figure 1.
- **Training cutoffs:** Boltz-2 2023-06-01; AlphaFold3 and OF3 p2 2021-09-30; OBv0 2025-06-30.

## Not included

The predicted structures (several GB per model, input and target) and the scripts that scored them are not in this repo; they live on MSKCC infrastructure. The processed tables in `data/` are what every figure is computed from, and `data/structure_ids.csv` lists every PDB entry involved so the experimental structures can be fetched from the PDB.

If you would like the actual predictions, get in touch with the authors and we will share them.

## License

Apache 2.0, see `LICENSE`.
