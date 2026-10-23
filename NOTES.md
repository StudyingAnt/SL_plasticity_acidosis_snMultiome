# Notes on the reorganized scripts

The scripts in this repository were reorganized from the analysis notebooks and HPC scripts used for the paper.

## Changes that do not affect results

- Hard-coded local and HPC paths were replaced by variables in `config/config.R` and `config/config.sh`.
- Notebook chunks were split into one script per analysis step, and every script reads its inputs from files written by the previous step.
- Exploratory analyses that are not part of the paper were removed. This includes ceramide generation axes, axis mean-expression and module scores, pairwise axis correlations, and alternative candidate-regulator plots.
- msVIPER was originally also run on two restricted gene universes: axis genes only, and KEGG sphingolipid genes. Only the full-transcriptome universe (`all_genes`) is used in the paper, so only that universe is kept. Each universe is computed independently, so dropping the others does not change the `all_genes` results.
- Plot-saving helpers and condition colors are defined once in `R/common.R`.
- `00-2_seurat_preprocessing/06_clustering_acdc.R` raises `future.globals.maxSize` to 10 GiB, the same value as in `05_transformations.R`. With recent versions of the `future` package, `FindNeighbors` on the merged object otherwise stops at the 500 MiB default. This setting is only a size limit.

## Changes in behavior

- `00-2_seurat_preprocessing/07_run_chromvar.R` now calls `set.seed(1234)` before `RunChromVAR()`. chromVAR samples background peaks at random, and the run for the paper had no fixed seed, so its deviation z-scores cannot be reproduced exactly. A rerun without a seed correlated with the original at r ≈ 0.97 (per-motif median r ≈ 0.97). In Figure 6, chromVAR is used only for the motif-support transparency of the regulator dot plots. The regulators shown, their order, NES and axis linkage do not depend on it.
- `02_aracne3/02_run_aracne3.sh` runs ARACNe3 inside the output directory and deletes `.ARACNe3_cached/` first. ARACNe3 caches its null model in the working directory and reuses it on later runs. A run that loads the cache skips the random draws used to build the null model, so its bootstrap subsamples differ from a fresh run with the same `--seed`. In our tests, a rerun with a cached null model shared only 43% of edges with the network used in the paper, whereas fresh runs contained all of its edges.
- `04_sphingolipid_axes/02_dominant_axis.R`: the `dominant4` label is now stored in the metacell object by matching metacell names. The original notebook assigned it by position after a `merge()` that re-orders rows. The dominant-axis UMAP and proportion panels are drawn directly from the merged table, so they are unaffected, and no script in this repository reads the stored `dominant4` column.

## Reproducibility check

The whole pipeline was rerun on a SLURM cluster (R 4.5.2, Bioconductor 3.22, Seurat 5.4.0, Signac 1.16.0, hdWGCNA 0.4.12), starting from the Cell Ranger ARC outputs, and compared with the objects and figures from the paper.

| Step | Result |
|---|---|
| 00-2 QC (01–02) | Nuclei per sample identical (9,318 / 12,774 / 9,821 / 9,442) |
| 00-2 clustering (06) | Same optimal parameters; `acdc_sct` and `acdc_peaks` labels identical for all 41,355 nuclei; `pca.sct` identical; SCT counts identical |
| 00-2 chromVAR (07) | r ≈ 0.97 with the original (random background peaks, see above) |
| 01 metacells | 5,387 metacells; CPM matrix byte-identical |
| 02 ARACNe3 | Regulator list identical (1,509 TFs). Contains 100.00% of the original edges (2,820,940 of 2,820,960) with identical mutual information, plus 0.27% additional edges |
| 03–05, Figure 6 | Every panel reproduced at the same size. Regulator dot plots show the same regulators on all four axes; two adjacent regulators swap order on the C1P and SM axes. Dominant-axis thresholds and diversity values differ from the paper only in the 3rd decimal place or later, and the ordering of conditions is unchanged for every measure |
