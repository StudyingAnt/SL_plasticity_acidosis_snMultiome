# Data availability

## Sequencing data (GEO)

Raw and processed single-nucleus multiome data for the four conditions are in GEO: **[GSE______](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE______)**.

| GEO sample | Library | Raw files (FASTQ) | Processed files |
|---|---|---|---|
| MCF7_Control_GEX / _ATAC | ct74 | SC60_1G (I1, I2, R1, R2) / SC60_1A (I1, R1, R2, R3) | `ct74_filtered_feature_bc_matrix.h5`, `ct74_per_barcode_metrics.csv` / `ct74_atac_fragments.tsv.gz` (+ `.tbi`) |
| MCF7_Acute_GEX / _ATAC | ct65 | SC60_2G / SC60_2A | `ct65_*` |
| MCF7_Chronic_GEX / _ATAC | aa65 | SC60_4G / SC60_4A | `aa65_*` |
| MCF7_Return_GEX / _ATAC | aa74 | SC60_3G / SC60_3A | `aa74_*` |

The processed files are Cell Ranger ARC 2.0.2 outputs, renamed with the sample prefix. To use them instead of rerunning [00-1_cellranger_arc](00-1_cellranger_arc), place them as `<DATASET_DIR>/<sample>/filtered_feature_bc_matrix.h5` and `<DATASET_DIR>/<sample>/atac_fragments.tsv.gz` (+ `.tbi`), and start from [00-2_seurat_preprocessing](00-2_seurat_preprocessing).

The series supplementary file `MCF7_snMultiome_cell_metadata.csv.gz` lists the 41,355 nuclei kept after QC. For each nucleus it gives the sample, condition, QC metrics and cluster labels (`acdc_sct`, `acdc_peaks`).

## Reference data

| Data | Source | Used in |
|---|---|---|
| GRCh38 reference `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` | 10x Genomics | 00-1 |
| Gene annotation EnsDb.Hsapiens.v86, genome BSgenome.Hsapiens.UCSC.hg38 | Bioconductor | 00-2 |
| JASPAR2024 CORE vertebrate motifs (*Homo sapiens*) | Bioconductor `JASPAR2024` | 00-2 step 07 |
| Human TFs v1.01 (`TF_names_v_1.01.txt`, 1,639 TFs) | https://humantfs.ccbr.utoronto.ca, downloaded and checksum-verified by `02_aracne3/00_download_human_tfs.sh` | 02 |

## Code

This repository: https://github.com/StudyingAnt/SL_plasticity_acidosis_snMultiome
