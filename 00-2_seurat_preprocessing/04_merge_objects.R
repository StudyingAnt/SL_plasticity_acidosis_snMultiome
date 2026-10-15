# =============================================================================
# 04. Merge the four per-sample objects.
# Output: <SEURAT_DIR>/all_qc_noMT_macs3.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))

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

cat(glue("[MERGE] Merging SeuratObjects across all samples...\n"))

merged <- merge(
  x = objs[[1]],
  y = objs[2:length(objs)],
  add.cell.ids = samples$sample,
  project = "all"
)

cat(glue(
  "[MERGE] Merge completed.\n",
  "  Assays: {paste(Assays(merged), collapse = ', ')}\n",
  "  Cells: {ncol(merged)}\n",
  "  Fragments: {length(Fragments(merged[['ATAC']]))}\n"
))

out_rds <- file.path(seurat_dir, "all_qc_noMT_macs3.rds")
cat(glue("[MERGE] Saving merged SeuratObject...\n"))
saveRDS(merged, out_rds)
cat(glue("[MERGE] Saved in: {out_rds}\n\n"))

cat(glue("[DONE] SeuratObject merge completed for all samples.\n"))
