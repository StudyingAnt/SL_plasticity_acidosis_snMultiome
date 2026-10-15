# Metacells

`01_build_metacells.R` aggregates transcriptionally similar nuclei into metacells with hdWGCNA. Low, sparse expression of sphingolipid metabolism genes makes per-nucleus estimates unreliable, which is why this step is needed.

| Parameter | Value |
|---|---|
| Input object | `all_qc_noMT_macs3_trf_clust.rds` (merged, 00-2 step 06) |
| Gene selection | `fraction`, 0.05 |
| Grouping | `condition` × `acdc_sct` |
| KNN space | `pca.sct` |
| Aggregation | sum of RNA raw counts |
| k | 25 |
| max_shared | 10 |
| min_cells | 30 |
| Seed | 1234 |

Genes with zero total count are removed, and the metacell count matrix is normalized to counts per million (CPM).

Outputs (`RESULTS_DIR/01_metacells`):

| File | Description |
|---|---|
| `all_metacell_hdWGCNA.rds` | Metacell Seurat object |
| `all_hdwgcna_object.rds` | Single-nucleus object with hdWGCNA data |
| `metacell_metadata.tsv` | Metacell metadata |
| `metacell_expression_counts.tsv` | Summed counts (genes × metacells) |
| `metacell_expression_CPM.tsv` | CPM, input to ARACNe3 and VIPER |
| `metacell_expression_log1pCPM.tsv` | log1p(CPM), not used downstream |
| `matrix_summary.txt` | Matrix size and library-size summary |
