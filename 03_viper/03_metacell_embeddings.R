# =============================================================================
# 03. Metacell embeddings:
#     (a) SCTransform on metacell RNA counts -> PCA (30) -> UMAP ("umap_sct_mc")
#     (b) Scaled VIPER activity -> PCA (30) -> UMAP ("umap_viper_mc")
#
# Input : <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER.rds
# Output: <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap.rds
#         <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(Seurat)
  library(SeuratObject)
})

source(file.path("R", "common.R"))

aracne_out_dir <- step_dir("03_viper")

in_rds <- file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER.rds")

if (!file.exists(in_rds)) {
  stop(paste("File not found:", in_rds))
}

metacell_obj <- readRDS(in_rds)

cat("Loaded:", in_rds, "\n")
cat("Class:", paste(class(metacell_obj), collapse = ", "), "\n")
cat("Cells:", ncol(metacell_obj), "\n")
cat("Assays:", paste(Assays(metacell_obj), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# (a) SCTransformation
# -----------------------------------------------------------------------------
options(future.globals.maxSize = 4 * 1024^3)

DefaultAssay(metacell_obj) <- "RNA"

# Recompute SCT on metacell RNA counts
metacell_obj <- SCTransform(
  metacell_obj,
  assay = "RNA",
  new.assay.name = "SCT_mc",
  verbose = FALSE
)

DefaultAssay(metacell_obj) <- "SCT_mc"

metacell_obj <- RunPCA(
  metacell_obj,
  assay = "SCT_mc",
  npcs = 30,
  reduction.name = "pca_sct_mc",
  verbose = FALSE
)

metacell_obj <- RunUMAP(
  metacell_obj,
  reduction = "pca_sct_mc",
  dims = 1:30,
  reduction.name = "umap_sct_mc",
  verbose = FALSE
)

saveRDS(
  metacell_obj,
  file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER_sct_umap.rds")
)

cat(glue("[metacell] Saved object with SCT-based metacell PCA/UMAP.\n"))

# -----------------------------------------------------------------------------
# (b) VIPER space
# -----------------------------------------------------------------------------
DefaultAssay(metacell_obj) <- "VIPER"

# Scale VIPER activity
metacell_obj <- ScaleData(
  metacell_obj,
  assay = "VIPER",
  features = rownames(metacell_obj[["VIPER"]]),
  verbose = FALSE
)

metacell_obj <- RunPCA(
  metacell_obj,
  assay = "VIPER",
  features = rownames(metacell_obj[["VIPER"]]),
  npcs = 30,
  reduction.name = "pca_viper_mc",
  verbose = FALSE
)

metacell_obj <- RunUMAP(
  metacell_obj,
  reduction = "pca_viper_mc",
  dims = 1:30,
  reduction.name = "umap_viper_mc",
  verbose = FALSE
)

saveRDS(
  metacell_obj,
  file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds")
)

cat(glue("[metacell] Saved object with VIPER PCA/UMAP.\n"))
