# Cell Ranger ARC

*Scripts for this step will be added.*

Raw FASTQ files were processed with Cell Ranger ARC (10x Genomics) against the GRCh38 reference. This produced, for each sample, the filtered feature-barcode matrix and the ATAC fragment file used in [00-2_seurat_preprocessing](../00-2_seurat_preprocessing).

Expected output layout (`DATASET_DIR` in `config/config.R`):

```
<DATASET_DIR>/
├── ct74/filtered_feature_bc_matrix.h5
├── ct74/atac_fragments.tsv.gz (+ .tbi)
├── ct65/...
├── aa65/...
└── aa74/...
```
