# =============================================================================
# 01. Build metacells with hdWGCNA and write the metacell expression matrices.
#
# Input : <SEURAT_DIR>/all_qc_noMT_macs3_trf_clust.rds (00-2_seurat_preprocessing)
# Output: <RESULTS_DIR>/01_metacells/
#           all_metacell_hdWGCNA.rds          metacell Seurat object
#           all_hdwgcna_object.rds            single-nucleus object with hdWGCNA slots
#           metacell_metadata.tsv
#           metacell_expression_counts.tsv    summed RNA counts
#           metacell_expression_CPM.tsv       CPM (input to ARACNe3 and VIPER)
#           metacell_expression_log1pCPM.tsv
#           matrix_summary.txt
# =============================================================================

source(file.path("R", "common.R"))

suppressPackageStartupMessages({
  library(glue)
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(Matrix)
  library(data.table)
})

set.seed(1234)

# =========================
# Paths
# =========================
merged_file <- file.path(SEURAT_DIR, "all_qc_noMT_macs3_trf_clust.rds")
out_dir <- step_dir("01_metacells")

# =========================
# Load original object
# =========================
cat(glue("[all] Loading merged SeuratObject from: {merged_file}\n"))
obj0 <- readRDS(merged_file)

cat(glue("[all] Original object loaded.\n"))
cat(glue("Cells: {ncol(obj0)}\n"))
cat(glue("Assays: {paste(names(obj0@assays), collapse = ', ')}\n"))
cat(glue("Reductions: {paste(Reductions(obj0), collapse = ', ')}\n"))

# =========================
# Create slim object for RNA/SCT only
# This avoids unnecessary ChromatinAssay/Seqinfo-related issues
# =========================
keep_assays <- intersect(c("RNA", "SCT"), names(obj0@assays))
keep_reductions <- intersect(c("pca.sct"), Reductions(obj0))

if (!("RNA" %in% keep_assays)) {
  stop("RNA assay not found in the object.")
}
if (!("SCT" %in% keep_assays)) {
  stop("SCT assay not found in the object.")
}
if (!("pca.sct" %in% keep_reductions)) {
  stop("Reduction 'pca.sct' not found in the object.")
}

obj <- DietSeurat(
  object = obj0,
  assays = keep_assays,
  dimreducs = keep_reductions,
  graphs = NULL,
  misc = TRUE
)

rm(obj0)
gc()

cat(glue("[all] Slim object created.\n"))
cat(glue("Cells: {ncol(obj)}\n"))
cat(glue("Assays: {paste(names(obj@assays), collapse = ', ')}\n"))
cat(glue("Reductions: {paste(Reductions(obj), collapse = ', ')}\n"))

# =========================
# Join RNA multilayers for Seurat v5
# =========================
obj <- JoinLayers(obj, assay = "RNA")

# =========================
# Metadata checks
# =========================
required_meta <- c("condition", "acdc_sct")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop(glue("Missing metadata columns: {paste(missing_meta, collapse = ', ')}"))
}

obj$condition <- as.character(obj$condition)
obj$acdc_sct <- as.character(obj$acdc_sct)

group_cols <- c("condition", "acdc_sct")
cat(glue("[all] group.by = {paste(group_cols, collapse = ', ')}\n"))

# =========================
# Setup hdWGCNA
# =========================
obj <- SetupForWGCNA(
  seurat_obj = obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "mcf7_metacell"
)

# =========================
# Construct metacells
# KNN is based on SCT-derived reduction
# Aggregation is based on RNA raw counts
# =========================
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = group_cols,
  reduction = "pca.sct",
  assay = "RNA",
  layer = "counts",
  mode = "sum",
  k = 25,
  max_shared = 10,
  min_cells = 30,
  ident.group = "acdc_sct"
)

# NormalizeMetacells is optional for hdWGCNA downstream
# It is not required for ARACNe3 CPM generation
# obj <- NormalizeMetacells(obj)

# =========================
# Extract metacell object
# =========================
metacell_obj <- GetMetacellObject(obj)
metacell_obj <- JoinLayers(metacell_obj, assay = "RNA")

cat(glue("[all] Metacell object created.\n"))
cat(glue("Metacells: {ncol(metacell_obj)}\n"))
cat(glue("Assays: {paste(names(metacell_obj@assays), collapse = ', ')}\n"))

# Save metacell object
metacell_rds <- file.path(out_dir, "all_metacell_hdWGCNA.rds")
saveRDS(metacell_obj, metacell_rds)
cat(glue("[all] Saved metacell object: {metacell_rds}\n"))

# =========================
# Save hdWGCNA-containing Seurat object
# =========================
hdwgcna_rds <- file.path(out_dir, "all_hdwgcna_object.rds")
saveRDS(obj, hdwgcna_rds)
cat(glue("[all] Saved hdWGCNA object: {hdwgcna_rds}\n"))

# =========================
# Extract RNA raw counts from metacells
# =========================
mc_counts <- LayerData(metacell_obj, assay = "RNA", layer = "counts")

if (is.null(mc_counts)) {
  stop("RNA counts layer not found in metacell object.")
}

cat(glue("[all] Extracted metacell RNA counts matrix: {nrow(mc_counts)} genes x {ncol(mc_counts)} metacells\n"))

# =========================
# Basic filtering
# Remove genes with zero total count across all metacells
# =========================
gene_sums <- Matrix::rowSums(mc_counts)
mc_counts <- mc_counts[gene_sums > 0, , drop = FALSE]

cat(glue("[all] Retained genes with nonzero total count: {nrow(mc_counts)}\n"))

# =========================
# Custom CPM normalization
# CPM = counts / column_sum * 1e6
# =========================
col_sums <- Matrix::colSums(mc_counts)

if (any(col_sums == 0)) {
  zero_cols <- names(col_sums)[col_sums == 0]
  stop(glue("Some metacells have zero total counts: {paste(zero_cols, collapse = ', ')}"))
}

cpm_mat <- t(t(mc_counts) / col_sums) * 1e6
cat(glue("[all] CPM matrix created: {nrow(cpm_mat)} genes x {ncol(cpm_mat)} metacells\n"))

# Optional log transform for sensitivity analysis
log1p_cpm_mat <- log1p(cpm_mat)

# =========================
# Save metacell metadata
# =========================
mc_meta <- metacell_obj@meta.data
mc_meta_file <- file.path(out_dir, "metacell_metadata.tsv")

mc_meta_dt <- data.table(
  metacell_id = rownames(mc_meta),
  mc_meta,
  check.names = FALSE
)

fwrite(
  mc_meta_dt,
  file = mc_meta_file,
  sep = "\t",
  quote = FALSE
)

cat(glue("[all] Saved metacell metadata: {mc_meta_file}\n"))

# =========================
# Save raw counts matrix
# =========================
counts_df <- data.frame(
  gene = rownames(mc_counts),
  as.matrix(mc_counts),
  check.names = FALSE
)

expr_counts_file <- file.path(out_dir, "metacell_expression_counts.tsv")
fwrite(
  counts_df,
  file = expr_counts_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(glue("[all] Saved raw counts matrix: {expr_counts_file}\n"))

# =========================
# Save ARACNe3 expression matrix: CPM
# ARACNe3 expects genes as rows and samples/metacells as columns
# =========================
cpm_df <- data.frame(
  gene = rownames(cpm_mat),
  as.matrix(cpm_mat),
  check.names = FALSE
)

expr_cpm_file <- file.path(out_dir, "metacell_expression_CPM.tsv")
fwrite(
  cpm_df,
  file = expr_cpm_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(glue("[all] Saved ARACNe3 CPM expression matrix: {expr_cpm_file}\n"))

# =========================
# Save optional log1p(CPM) matrix
# =========================
log1p_cpm_df <- data.frame(
  gene = rownames(log1p_cpm_mat),
  as.matrix(log1p_cpm_mat),
  check.names = FALSE
)

expr_log1p_cpm_file <- file.path(out_dir, "metacell_expression_log1pCPM.tsv")
fwrite(
  log1p_cpm_df,
  file = expr_log1p_cpm_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat(glue("[all] Saved optional log1p(CPM) matrix: {expr_log1p_cpm_file}\n"))

# The ARACNe3 regulator list is prepared in 02_aracne3/01_prepare_regulators.R

# =========================
# Quick QC summaries
# =========================
qc_file <- file.path(out_dir, "matrix_summary.txt")
writeLines(
  c(
    paste0("metacells=", ncol(mc_counts)),
    paste0("genes_nonzero=", nrow(mc_counts)),
    paste0("min_metacell_total_counts=", min(col_sums)),
    paste0("median_metacell_total_counts=", median(col_sums)),
    paste0("max_metacell_total_counts=", max(col_sums))
  ),
  con = qc_file
)

cat(glue("[all] Saved matrix summary: {qc_file}\n"))
cat("[all] Done.\n")
