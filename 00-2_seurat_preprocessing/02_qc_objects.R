# =============================================================================
# 02. Filter nuclei by QC thresholds (config/qc_parameters.json) and remove
#     mitochondrial genes from the RNA assay.
# Output: <SEURAT_DIR>/<sample>_qc_noMT.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))

# QC parameters json
qc_json_path <- file.path("config", "qc_parameters.json")

if (!file_exists(qc_json_path)) {
  stop(glue("QC parameter file not found: {qc_json_path}"))
}

qc <- fromJSON(qc_json_path)

nCount_ATAC_min <- qc$nCount_ATAC$min
nCount_ATAC_max <- qc$nCount_ATAC$max
nucleosome_signal_max <- qc$nucleosome_signal$max
TSS_enrichment_min <- qc$TSS_enrichment$min
nCount_RNA_min <- qc$nCount_RNA$min
nCount_RNA_max <- qc$nCount_RNA$max
percent_mt_max <- qc$percent_mt$max

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]

  # Load raw SeuratObjects
  obj_file <- file.path(seurat_dir, paste0(smp, "_raw.rds"))
  cat(glue("[{smp}] Loading SeuratObject from: {obj_file}...\n"))
  obj <- readRDS(obj_file)

  cat(glue(
    "[{smp}] Filtering cells using QC parameters:\n",
    "  nCount_ATAC > {nCount_ATAC_min} & nCount_ATAC < {nCount_ATAC_max}\n",
    "  nucleosome_signal < {nucleosome_signal_max}\n",
    "  TSS.enrichment > {TSS_enrichment_min}\n",
    "  nCount_RNA > {nCount_RNA_min} & nCount_RNA < {nCount_RNA_max}\n",
    "  percent.mt < {percent_mt_max}\n"
  ))

  n_cells_before <- ncol(obj)

  qc_obj <- subset(
    x = obj,
    subset =
      nCount_ATAC < nCount_ATAC_max &
      nCount_ATAC > nCount_ATAC_min &
      nucleosome_signal < nucleosome_signal_max &
      TSS.enrichment > TSS_enrichment_min &
      nCount_RNA < nCount_RNA_max &
      nCount_RNA > nCount_RNA_min &
      percent.mt < percent_mt_max
  )

  n_cells_after <- ncol(qc_obj)
  n_cells_removed <- n_cells_before - n_cells_after
  pct_kept <- 100 * n_cells_after / n_cells_before
  pct_removed <- 100 * n_cells_removed / n_cells_before

  cat(glue(
    "[{smp}] QC filtering summary:\n",
    "  Before: {n_cells_before} cells\n",
    "  After: {n_cells_after} cells\n",
    "  Removed: {n_cells_removed} cells ({sprintf('%.2f', pct_removed)}%)\n",
    "  Kept: {sprintf('%.2f', pct_kept)}%\n"
  ))

  cat(glue("[{smp}] Removing mitochondrial genes from RNA assay...\n"))
  rna_counts <- LayerData(qc_obj, assay = "RNA", layer = "counts")
  mt_idx <- grep("^MT-", rownames(rna_counts))

  if (length(mt_idx) == 0) {
    cat(glue("[{smp}] No mitochondrial genes matched '^MT-'. Keeping all RNA features.\n"))
    rna_counts_noMT <- rna_counts
  } else {
    rna_counts_noMT <- rna_counts[-mt_idx, ]
  }

  proj_name <- paste0(smp, "_", cond)

  qc_noMT_obj <- CreateSeuratObject(
    counts = rna_counts_noMT,
    assay = "RNA",
    project = proj_name,
    meta.data = qc_obj@meta.data
  )

  qc_noMT_obj[["ATAC"]] <- qc_obj[["ATAC"]]

  cat(glue(
    "[{smp}] noMT object created.\n",
    "  Condition: {cond}\n",
    "  Cells: {ncol(qc_noMT_obj)}\n",
    "  RNA features: {nrow(qc_noMT_obj[['RNA']])}\n",
    "  ATAC features: {nrow(qc_noMT_obj[['ATAC']])}\n",
    "  Fragments: {length(Fragments(qc_noMT_obj[['ATAC']]))}\n"
  ))

  qc_noMT_file <- file.path(seurat_dir, paste0(smp, "_qc_noMT.rds"))
  cat(glue("[{smp}] Saving noMT SeuratObject as RDS...\n"))
  saveRDS(qc_noMT_obj, qc_noMT_file)
  cat(glue("[{smp}] Saved in: {qc_noMT_file}\n"))
}

cat(glue("[DONE] QC filtering and noMT SeuratObject creation completed for all samples.\n"))
