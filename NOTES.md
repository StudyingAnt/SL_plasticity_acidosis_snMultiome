# Notes on the reorganized scripts

The scripts in this repository were reorganized from the analysis notebooks and HPC scripts used for the paper.

## Changes that do not affect results

- Hard-coded local and HPC paths were replaced by variables in `config/config.R` and `config/config.sh`.
- Notebook chunks were split into one script per analysis step, and every script reads its inputs from files written by the previous step.
- Exploratory analyses that are not part of the paper were removed. This includes ceramide generation axes, axis mean-expression and module scores, pairwise axis correlations, and alternative candidate-regulator plots.
- msVIPER was originally also run on two restricted gene universes: axis genes only, and KEGG sphingolipid genes. Only the full-transcriptome universe (`all_genes`) is used in the paper, so only that universe is kept. Each universe is computed independently, so dropping the others does not change the `all_genes` results.
- Plot-saving helpers and condition colors are defined once in `R/common.R`.

## Changes in behavior

- `04_sphingolipid_axes/02_dominant_axis.R`: the `dominant4` label is now stored in the metacell object by matching metacell names. The original notebook assigned it by position after a `merge()` that re-orders rows. Figures 6c and 6d are drawn directly from the merged table, so they are unaffected, and no script in this repository reads the stored `dominant4` column.
