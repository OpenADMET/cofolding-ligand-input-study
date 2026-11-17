# Data

Only what the figures in the post use, already processed. The predicted structures and the raw per-run score files live on MSKCC infrastructure and are not part of this repo.

| file | used by | what it is |
|---|---|---|
| `structure_ids.csv` | all | **every PDB entry behind the figures**: `set`, `pdb_id`, release date, and the figures it appears in. Regenerate with `scripts/make_structure_ids.py`. |
| `pxr_cyp3a4_pose_rmsds.csv.gz` | Figs 1, 3 (AlphaFold3), 4 | one row per scored pose for the 49 PXR and 85 CYP3A4 structures: `model`, `target` (pxr / cyp3a4), `input` (smiles / ccd / ccd_genericnames), `seed`, `structure`, `sample`, `ligand_rmsd`, `status`. Boltz-2, AlphaFold3, OF3 (p2) and OBv0, five seeds x five samples, no templates. |
| `pdb_pxr_cyp3a4_release_dates.csv` | Fig 4 | release date of each of those PDB entries, for the training-cutoff split |
| `names_effect_per_structure.csv` | Fig 3 | per structure, whether the best of 5 samples was within 2 Å at seed 42 under each step of the atom-name experiment (SMILES, generic names, deposited names), for Boltz-2 / OF3 / OBv0 |
| `random_complexes_before_cutoff.csv` | Fig 5 | poses for the 50 random PDB complexes released 2017-2021: `model`, `pdb_id`, `mode` (smiles / ccd), `sample_idx`, `BiSyRMSD` |
| `random_complexes_after_cutoff.csv` | Fig 5 | same, for the 50 complexes released after 2025-06-30 |
| `random_complexes_info.csv` | `structure_ids.csv` | release date and ligand CCD code of each of those 100 complexes |

Notes
- In the two `random_complexes_*` files the OpenFold model is labelled `of3`; the post shows it as OBv0 (see Fig 5).
- In `names_effect_per_structure.csv` the model `OF3 p2` is shown as OF3.
- The three ligand inputs are `smiles`, `ccd` (the CCD code), and a generic-atom-names arm (`ccd_genericnames` for AlphaFold3 via a user CCD; the stored-conformer / isomeric-SMILES equivalents for Boltz-2 and the OpenFold models, as described in the post).
