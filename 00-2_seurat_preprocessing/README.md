# Seurat / Signac preprocessing

These scripts build the per-sample and merged Seurat objects from the Cell Ranger ARC output. They were run on an HPC cluster (Linux). `00_setup.R` is sourced by every script and loads packages, the sample table and the hg38 annotation (EnsDb.Hsapiens.v86, UCSC style).

| Script | Description | Output (`SEURAT_DIR`) |
|---|---|---|
| `01_make_objects.R` | RNA assay + ATAC ChromatinAssay per sample; nucleosome signal, TSS enrichment, % mitochondrial reads | `<sample>_raw.rds` |
| `02_qc_objects.R` | QC filtering, removal of mitochondrial genes | `<sample>_qc_noMT.rds` |
| `03_call_peaks.R` | MACS3 peaks per sample, standard chromosomes only, hg38 blacklist removed, union peak set, unified `peaks` assay | `<sample>_qc_noMT_macs3.rds` |
| `04_merge_objects.R` | Merge the four samples | `all_qc_noMT_macs3.rds` |
| `05_transformations.R` | SCTransform → PCA (50) → UMAP; TF-IDF → SVD (LSI, dims 2:51) → UMAP; WNN | `<name>_qc_noMT_macs3_trf.rds` |
| `06_clustering_acdc.R` | ACDC-style grid-search Leiden clustering on `pca.sct` and `peaks.lsi` | `<name>_qc_noMT_macs3_trf_clust.rds` |
| `07_run_chromvar.R` | JASPAR2024 motifs and chromVAR on the merged object | `all_qc_noMT_macs3_trf_chrvar.rds` |

Downstream steps use `all_qc_noMT_macs3_trf_clust.rds` and `all_qc_noMT_macs3_trf_chrvar.rds`.

## QC thresholds

`config/qc_parameters.json`

| Metric | Kept |
|---|---|
| ATAC fragment count | 1,000 – 100,000 |
| Nucleosome signal | < 2 |
| TSS enrichment | > 1 |
| RNA UMI count | 1,000 – 35,000 |
| Mitochondrial reads | < 20% |

## Clustering grid

`06_clustering_acdc.R` takes one object per call (`-i`):

| Parameter | Value |
|---|---|
| k (nearest neighbors) | 11, 21, …, 101 |
| Resolution | 0.1, 0.3, …, 1.9 |
| Algorithm | Leiden (`algorithm = 4`) |
| Criterion | Mean of per-cluster mean silhouette widths (Euclidean, PCA or LSI space) |

The selected clusters are stored as `acdc_sct` (RNA) and `acdc_peaks` (ATAC).

```bash
Rscript 00-2_seurat_preprocessing/06_clustering_acdc.R -i <SEURAT_DIR>/all_qc_noMT_macs3_trf.rds
```

## chromVAR

Motif position frequency matrices come from JASPAR2024 (CORE collection, vertebrates, *Homo sapiens*, latest versions only). Motifs are matched to peaks with BSgenome.Hsapiens.UCSC.hg38.
