# Extended Data: ATAC overview

| Script | Description | Panels |
|---|---|---|
| `01_plot_atac_overview.R` | UMAPs of all nuclei by condition (RNA, ATAC, WNN) and chromVAR motif deviations of the Figure 6f–i candidate regulators | a–d |

## Panels

| Panel | Content |
|---|---|
| a | UMAP from RNA (SCT PCA, `rna.sct.umap`) |
| b | UMAP from ATAC (peaks LSI, `peaks.umap`) |
| c | UMAP from WNN (RNA + ATAC, `mm.sct.peaks.umap`) |
| d | Mean per-nucleus chromVAR deviation (z) per condition for the JASPAR2024 motifs of the candidate regulators |

Panel d uses the same regulators as the Figure 6f–i dot plots (msVIPER padj < 0.05, top 20 per axis by mean −log10 padj) and the same TF → motif name matching as `05_axis_regulators/02_regulon_linkage_motif_support.R`. Regulators without a JASPAR motif are not shown.

## Inputs

- `<SEURAT_DIR>/all_qc_noMT_macs3_trf_chrvar.rds` (00-2 step 07)
- `<RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/viz_allgenes_regulon_motifalpha_padj005/plot_dt_allgenes_regulon_motifalpha_padj005.tsv` (step 5)

## Outputs

Written to `RESULTS_DIR/06_extended_data_atac`:

- `umap_{rna,atac,wnn}_by_condition.{svg,png}`, `umap_legend.{svg,png}`
- `chromvar_candidate_regulators_heatmap.{svg,png}`
- `chromvar_candidate_regulators_mean_by_condition.tsv` (values plotted in panel d)
