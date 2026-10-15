# =============================================================================
# 05. RNA: SCTransform -> PCA (50) -> UMAP; ATAC peaks: TF-IDF -> SVD (LSI,
#     dims 2:51) -> UMAP; WNN integration. Run on each sample and on the merged
#     object.
# Output: <SEURAT_DIR>/<name>_qc_noMT_macs3_trf.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))
suppressMessages({
  library(future)
})

options(future.globals.maxSize = 10 * 1024^3)

objs <- list()

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]

  obj_file <- file.path(seurat_dir, paste0(smp, "_qc_noMT_macs3.rds"))
  cat(glue("[{smp}] Loading SeuratObject from: {obj_file}\n"))
  obj <- readRDS(obj_file)

  cat(glue(
    "[{smp}] SeuratObject loaded.\n",
    "  Condition: {cond}\n",
    "  Assays: {paste(Assays(obj), collapse = ', ')}\n",
    "  Cells: {ncol(obj)}\n",
    "  Fragments: {length(Fragments(obj[['ATAC']]))}\n"
  ))

  objs[[smp]] <- obj
}

merged_file <- file.path(seurat_dir, "all_qc_noMT_macs3.rds")
cat(glue("[all] Loading merged SeuratObject from: {merged_file}\n"))
merged <- readRDS(merged_file)

cat(glue(
  "[all] SeuratObject loaded.\n",
  "  Condition: merged\n",
  "  Assays: {paste(Assays(merged), collapse = ', ')}\n",
  "  Cells: {ncol(merged)}\n",
  "  Fragments: {length(Fragments(merged[['ATAC']]))}\n"
))

objs[["all"]] <- merged

for (obj_name in names(objs)) {
  if (obj_name %in% samples$sample) {
    cond <- samples$condition[match(obj_name, samples$sample)]
  } else {
    cond <- "merged"
  }

  cat(glue("[{obj_name}] Running SCTransform on RNA assay...\n"))

  obj <- objs[[obj_name]]

  obj <- SCTransform(
    object = obj,
    assay = "RNA"
  )

  obj <- RunPCA(
    object = obj,
    assay = "SCT",
    npcs = 50,
    reduction.name = "pca.sct",
    reduction.key = "RNASCTPCA"
  )

  obj <- FindNeighbors(
    object = obj,
    reduction = "pca.sct",
    dims = 1:50,
    graph.name = c("SCT.nn", "SCT.snn")
  )

  obj <- RunUMAP(
    object = obj,
    reduction = "pca.sct",
    dims = 1:50,
    metric = "cosine",
    reduction.name = "rna.sct.umap",
    reduction.key = "RNASCTUMAP"
  )

  cat(glue(
    "[{obj_name}] RNA transformation completed.\n",
    "  Condition: {cond}\n",
    "  Assays: {paste(Assays(obj), collapse = ', ')}\n",
    "  Cells: {ncol(obj)}\n"
  ))

  objs[[obj_name]] <- obj
}

for (obj_name in names(objs)) {
  if (obj_name %in% samples$sample) {
    cond <- samples$condition[match(obj_name, samples$sample)]
  } else {
    cond <- "merged"
  }

  cat(glue("[{obj_name}] Running LSI transformation on peaks assay...\n"))

  obj <- objs[[obj_name]]

  obj <- FindTopFeatures(
    object = obj,
    assay = "peaks",
    min.cutoff = 5
  )

  obj <- RunTFIDF(
    object = obj,
    assay = "peaks"
  )

  obj <- RunSVD(
    object = obj,
    assay = "peaks",
    n = 51,
    reduction.name = "peaks.lsi",
    reduction.key = "PEAKSLSI"
  )

  obj <- FindNeighbors(
    object = obj,
    reduction = "peaks.lsi",
    dims = 2:51,
    graph.name = c("peaks.nn", "peaks.snn")
  )

  obj <- RunUMAP(
    object = obj,
    reduction = "peaks.lsi",
    dims = 2:51,
    metric = "cosine",
    reduction.name = "peaks.umap",
    reduction.key = "PEAKSUMAP"
  )

  cat(glue(
    "[{obj_name}] Peaks transformation completed.\n",
    "  Condition: {cond}\n",
    "  peaks features: {nrow(obj[['peaks']])}\n",
    "  Cells: {ncol(obj)}\n"
  ))

  objs[[obj_name]] <- obj
}

for (obj_name in names(objs)) {
  if (obj_name %in% samples$sample) {
    cond <- samples$condition[match(obj_name, samples$sample)]
  } else {
    cond <- "merged"
  }

  cat(glue("[{obj_name}] Running multimodal integration...\n"))

  obj <- objs[[obj_name]]

  obj <- FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("pca.sct", "peaks.lsi"),
    dims.list = list(1:50, 2:51),
    knn.graph.name = "peaks.sct.wknn",
    snn.graph.name = "peaks.sct.wsnn",
    weighted.nn.name = "peaks.sct.weighted.nn",
    modality.weight.name = c("SCT.weight", "peaks.sct.weight")
  )

  obj <- RunUMAP(
    object = obj,
    metric = "cosine",
    nn.name = "peaks.sct.weighted.nn",
    reduction.name = "mm.sct.peaks.umap"
  )

  cat(glue(
    "[{obj_name}] Multimodal integration completed.\n",
    "  Condition: {cond}\n",
    "  Reductions: {paste(Reductions(obj), collapse = ', ')}\n",
    "  Cells: {ncol(obj)}\n"
  ))

  objs[[obj_name]] <- obj
}

for (obj_name in names(objs)) {
  obj <- objs[[obj_name]]

  out_rds <- file.path(seurat_dir, paste0(obj_name, "_qc_noMT_macs3_trf.rds"))
  cat(glue("[{obj_name}] Saving transformed SeuratObject...\n"))
  saveRDS(obj, out_rds)
  cat(glue("[{obj_name}] Saved in: {out_rds}\n\n"))
}

cat(glue(
  "[DONE] RNA transformation, peaks LSI transformation, multimodal integration, and SeuratObject saving completed for all samples and the merged object.\n"
))
