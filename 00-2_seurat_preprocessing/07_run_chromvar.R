# =============================================================================
# 07. Add JASPAR2024 CORE vertebrate (human) motifs to the peaks assay and run
#     chromVAR on the merged object.
# Output: <SEURAT_DIR>/all_qc_noMT_macs3_trf_chrvar.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))
suppressMessages({
  library(TFBSTools)
  library(JASPAR2024)
})

merged_file <- file.path(seurat_dir, "all_qc_noMT_macs3_trf_clust.rds")
cat(glue("[all] Loading merged SeuratObject from: {merged_file}\n"))
obj <- readRDS(merged_file)

cat(glue(
  "[all] SeuratObject loaded.\n",
  "  Condition: merged\n",
  "  Assays: {paste(Assays(obj), collapse = ', ')}\n",
  "  Cells: {ncol(obj)}\n",
  "  Fragments: {length(Fragments(obj[['ATAC']]))}\n"
))

cat(glue("[JASPAR2024] Loading motif PFMs from JASPAR2024...\n"))

jaspar_db <- JASPAR2024::JASPAR2024()
jaspar_con <- RSQLite::dbConnect(
  RSQLite::SQLite(),
  JASPAR2024::db(jaspar_db)
)

pfm <- TFBSTools::getMatrixSet(
  x = jaspar_con,
  opts = list(
    collection = "CORE",
    tax_group = "vertebrates",
    species = 9606,
    all_versions = FALSE
  )
)


#jaspar_db <- JASPAR2024::JASPAR2024()
#pfm <- TFBSTools::getMatrixSet(
#  x = jaspar_db,
#  opts = list(
#    collection = "CORE",
#    tax_group = "vertebrates",
#    species = 9606,
#    all_versions = FALSE
#  )
#)

cat(glue(
  "[JASPAR2024] Motif PFMs loaded.\n",
  "  Number of PFMs: {length(pfm)}\n"
))

cat(glue("[all] Adding motif information to peaks assay...\n"))

obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "peaks"
)

cat(glue(
  "[all] Motif information added.\n",
  "  Condition: merged\n",
  "  peaks features: {nrow(obj[['peaks']])}\n"
))

cat(glue("[all] Running chromVAR on peaks assay...\n"))

# chromVAR samples background peaks at random; fix the seed so reruns give
# identical deviations (the run for the paper had no seed, see NOTES.md)
set.seed(1234)

obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "peaks"
)

cat(glue(
  "[all] chromVAR completed.\n",
  "  Condition: merged\n",
  "  Assays: {paste(Assays(obj), collapse = ', ')}\n"
))

out_rds <- file.path(seurat_dir, "all_qc_noMT_macs3_trf_chrvar.rds")
cat(glue("[all] Saving SeuratObject with motifs and chromVAR...\n"))
saveRDS(obj, out_rds)
cat(glue("[all] Saved in: {out_rds}\n\n"))

RSQLite::dbDisconnect(jaspar_con)

cat(glue(
  "[DONE] JASPAR2024 motif loading, motif attachment to peaks assay, chromVAR, and SeuratObject saving completed for the merged object.\n"
))
