# =============================================================================
# 04. Attach metacell-level chromVAR deviations: average the per-nucleus
#     chromVAR z-scores over the member nuclei of each metacell
#     (hdWGCNA "cells_merged"), and store JASPAR motif names in @misc.
#
# Input : <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds
#         <SEURAT_DIR>/all_qc_noMT_macs3_trf_chrvar.rds (00-2_seurat_preprocessing/07)
# Output: <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap_chvar.rds
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(Matrix)
})

source(file.path("R", "common.R"))

aracne_out_dir <- step_dir("03_viper")

metacell_obj <- readRDS(file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds"))

merged_file <- file.path(SEURAT_DIR, "all_qc_noMT_macs3_trf_chrvar.rds")

cat(glue("[all] Loading merged SeuratObject from: {merged_file}
"))
obj0 <- readRDS(merged_file)

add_metacell_chromvar_assay <- function(
  obj0,
  metacell_obj,
  chromvar_assay = "chromvar",
  chromvar_layer = "data",
  membership_col = "cells_merged",
  output_assay = "chromvar_mc",
  verbose = TRUE
) {
  if (!(chromvar_assay %in% names(obj0@assays))) {
    stop(paste0("Assay not found in obj0: ", chromvar_assay))
  }

  if (!(membership_col %in% colnames(metacell_obj@meta.data))) {
    stop(paste0("Column not found in metacell_obj@meta.data: ", membership_col))
  }

  chrom_mat <- GetAssayData(
    object = obj0,
    assay = chromvar_assay,
    layer = chromvar_layer
  )

  if (is.null(chrom_mat)) {
    stop("chromVAR matrix is NULL.")
  }

  if (verbose) {
    cat("[chromVAR] Loaded chromVAR matrix from obj0.\n")
    cat("Features:", nrow(chrom_mat), "\n")
    cat("Cells:", ncol(chrom_mat), "\n")
  }

  mc_ids <- colnames(metacell_obj)
  membership_list <- metacell_obj@meta.data[[membership_col]]

  if (length(membership_list) != length(mc_ids)) {
    stop("Length of membership column does not match number of metacells.")
  }

  mc_chrom_list <- vector("list", length = length(mc_ids))
  names(mc_chrom_list) <- mc_ids

  n_missing_total <- 0L

  for (i in seq_along(mc_ids)) {
    mc_id <- mc_ids[i]
    member_cells <- membership_list[[i]]

    if (is.null(member_cells) || length(member_cells) == 0) {
      stop(paste0("Metacell has no member cells: ", mc_id))
    }

    member_cells <- as.character(member_cells)
    keep_cells <- intersect(member_cells, colnames(chrom_mat))
    n_missing <- length(setdiff(member_cells, keep_cells))
    n_missing_total <- n_missing_total + n_missing

    if (length(keep_cells) == 0) {
      stop(paste0("No member cells found in chromVAR matrix for metacell: ", mc_id))
    }

    mc_vec <- Matrix::rowMeans(chrom_mat[, keep_cells, drop = FALSE])
    mc_chrom_list[[i]] <- mc_vec

    if (verbose && (i %% 100 == 0 || i == length(mc_ids))) {
      cat("[chromVAR] Processed", i, "/", length(mc_ids), "metacells\n")
    }
  }

  mc_chrom_mat <- do.call(cbind, mc_chrom_list)
  colnames(mc_chrom_mat) <- mc_ids

  if (!all(colnames(mc_chrom_mat) == colnames(metacell_obj))) {
    stop("Final metacell chromVAR matrix columns do not match metacell_obj columns.")
  }

  metacell_obj[[output_assay]] <- CreateAssayObject(data = mc_chrom_mat)

  if (verbose) {
    cat("[chromVAR] Added assay to metacell_obj:", output_assay, "\n")
    cat("[chromVAR] Matrix size:", nrow(mc_chrom_mat), "features x", ncol(mc_chrom_mat), "metacells\n")
    cat("[chromVAR] Total missing member-cell matches skipped:", n_missing_total, "\n")
  }

  metacell_obj
}


metacell_obj@meta.data$cells_merged <- lapply(
  metacell_obj@meta.data$cells_merged,
  function(x) trimws(strsplit(as.character(x), split = ",", fixed = TRUE)[[1]])
)

metacell_obj <- add_metacell_chromvar_assay(
  obj0 = obj0,
  metacell_obj = metacell_obj,
  chromvar_assay = "chromvar",
  chromvar_layer = "data",
  membership_col = "cells_merged",
  output_assay = "chromvar_mc",
  verbose = TRUE
)

# Extract motif.names from obj0 peaks motifs
motif_names <- obj0@assays[["peaks"]]@motifs@motif.names
if (is.null(motif_names) || length(motif_names) == 0) {
  stop("obj0@assays[['peaks']]@motifs@motif.names is NULL/empty.")
}

# Store into metacell_obj@misc as 'motif.names' (exact name)
metacell_obj@misc$motif.names <- motif_names

cat("Stored metacell_obj@misc$motif.names\n")
cat("Length:", length(metacell_obj@misc$motif.names), "\n")
cat("Example motif IDs:", paste(head(names(metacell_obj@misc$motif.names), 3), collapse = ", "), "\n")

saveRDS(
  metacell_obj,
  file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap_chvar.rds")
)

cat(glue("[metacell] Saved object with metacell chromVAR assay.\n"))
