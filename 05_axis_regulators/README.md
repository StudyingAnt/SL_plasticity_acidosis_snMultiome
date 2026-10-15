# Candidate regulators of each sphingolipid axis

| Script | Description | Figure |
|---|---|---|
| `01_msviper_topbottom.R` | msVIPER for each condition × axis contrast | – |
| `02_regulon_linkage_motif_support.R` | Adds axis gene linkage and chromVAR motif support | – |
| `03_plot_regulator_dotplots.R` | Dot plots of candidate regulators | 6g–j |

## msVIPER

For each condition and axis:

1. Metacells are ranked by axis pressure, and the top 20% and bottom 20% are selected.
2. A gene signature is computed on the full-transcriptome CPM matrix with row-wise t-tests (`genefilter::rowttests`). Two-sided p-values are converted to signed z-scores.
3. msVIPER is run with `minsize = 25` and `ges.filter = TRUE`.

## Evidence integration

| Evidence | Definition |
|---|---|
| msVIPER | NES and Benjamini–Hochberg adjusted p-value |
| Axis gene linkage | Σ \|tfmode × likelihood\| from the regulator to the axis genes in the regulon |
| Motif support | Wilcoxon test of metacell chromVAR deviations (top vs bottom 20%) for the regulator's JASPAR2024 motif(s); best BH-adjusted p-value |

## Dot plots

The dot plots include regulators with msVIPER padj < 0.05. For each axis, the top 20 regulators are chosen by mean −log10(padj) across conditions.

| Encoding | Value |
|---|---|
| Color | NES |
| Size | Axis gene linkage |
| Transparency | −log10(motif padj) |

Outputs are written to `RESULTS_DIR/05_axis_regulators/msviper_topbottom20_by_axis`.
