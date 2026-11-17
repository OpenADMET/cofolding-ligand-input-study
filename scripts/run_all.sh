#!/usr/bin/env bash
# Re-run every figure notebook in place; the PNGs land in figures/.
set -euo pipefail
cd "$(dirname "$0")/../notebooks"
for nb in fig1_smiles_vs_ccd fig2_input_schematic fig3_atom_names fig4_cutoff_split fig5_random_complexes; do
  echo "== $nb"
  jupyter nbconvert --to notebook --execute --inplace "$nb.ipynb"
done
