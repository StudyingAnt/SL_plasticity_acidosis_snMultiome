# =============================================================================
# 01. Regulatory pressure score for each ceramide clearance axis.
#
#   P_gs = sum_r(w_rg * A_rs) / sum_r(|w_rg|),  w_rg = tfmode_rg * likelihood_rg
#   M_s  = mean of P_gs over the axis gene set
#   (A = VIPER activity, w = ARACNe3/VIPER regulon weights)
#
# Input : <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap_chvar.rds
#         <RESULTS_DIR>/03_viper/mcf7_metacell_regulon_viper.rds
# Output: <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores.rds
#         (meta.data columns Pressure_s1p, Pressure_c1p, Pressure_sm, Pressure_gsl)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

source(file.path("R", "common.R"))

out_dir <- step_dir("04_sphingolipid_axes")

metacell_obj <- readRDS(file.path(DIR_VIPER, "all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap_chvar.rds"))
regulon_obj <- readRDS(file.path(DIR_VIPER, "mcf7_metacell_regulon_viper.rds"))

metacell_obj$condition <- factor(
  metacell_obj$condition,
  levels = CONDITION_LEVELS
)

compute_target_pressure <- function(
  regulon_obj,
  viper_mat,
  targets,
  regulators = NULL,
  layer_name = NULL,
  eps = 1e-12,
  sparse_W = TRUE,
  return_parts = FALSE
) {
  # viper_mat can be a matrix (regulator x cell) or a Seurat object if layer_name is provided
  if (!is.null(layer_name)) {
    # layer_name is only for Seurat v5 workflow; pass viper_mat as Seurat object in that case
    viper_mat <- SeuratObject::GetAssayData(viper_mat, assay = "VIPER", layer = layer_name)
  }

  A <- as.matrix(viper_mat)
  if (is.null(rownames(A))) stop("viper_mat must have rownames (regulators).")
  if (is.null(colnames(A))) stop("viper_mat must have colnames (cells/samples).")

  if (is.null(regulators)) {
    regulators <- intersect(rownames(A), names(regulon_obj))
  } else {
    regulators <- intersect(regulators, intersect(rownames(A), names(regulon_obj)))
  }

  if (length(regulators) == 0) {
    stop("No overlapping regulators among viper_mat rownames, regulators, and regulon_obj names.")
  }

  targets <- unique(targets)
  if (length(targets) == 0) stop("targets is empty.")

  # Build W (regulator x target) where W[r,g] = tfmode_rg * likelihood_rg
  W <- matrix(
    0,
    nrow = length(regulators),
    ncol = length(targets),
    dimnames = list(regulators, targets)
  )

  for (reg in regulators) {
    entry <- regulon_obj[[reg]]
    tfmode_vec <- entry$tfmode
    likelihood_vec <- entry$likelihood

    if (is.null(tfmode_vec) || is.null(likelihood_vec)) next

    idx <- match(targets, names(tfmode_vec))
    ok <- which(!is.na(idx))
    if (length(ok) == 0) next

    # Use tfmode names to locate targets; use same indices for likelihood (aligned vectors)
    tfv <- unname(tfmode_vec[idx[ok]])
    lkv <- unname(likelihood_vec[idx[ok]])
    W[reg, ok] <- tfv * lkv
  }

  # Restrict A to selected regulators and align ordering with W
  A2 <- A[regulators, , drop = FALSE]

  # Numerator: target x cell = t(W) %*% A
  if (sparse_W) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required when sparse_W = TRUE.")
    }
    Wm <- Matrix::Matrix(W, sparse = TRUE)
    numerator <- as.matrix(Matrix::t(Wm) %*% A2)
  } else {
    numerator <- t(W) %*% A2
  }

  # Denominator: target vector = sum_r |w_rg|
  denominator <- colSums(abs(W))
  denominator[denominator < eps] <- NA_real_

  # Pressure: divide each target row by its denominator
  pressure <- sweep(numerator, 1, denominator, "/")

  # Keep dimnames
  rownames(pressure) <- targets
  colnames(pressure) <- colnames(A2)
  rownames(numerator) <- targets
  colnames(numerator) <- colnames(A2)
  names(denominator) <- targets

  if (return_parts) {
    return(list(
      pressure = pressure,
      numerator = numerator,
      denominator = denominator,
      W = W,
      regulators_used = regulators,
      targets_used = targets
    ))
  }

  return(pressure)
}

compute_axis_score <- function(P, genes, na.rm = TRUE, warn_missing = TRUE) {
  # Check input matrix
  P <- as.matrix(P)
  if (is.null(rownames(P))) {
    stop("P must have rownames corresponding to target genes.")
  }

  genes <- unique(genes)
  if (length(genes) == 0) {
    stop("genes is empty.")
  }

  present_genes <- intersect(genes, rownames(P))
  missing_genes <- setdiff(genes, rownames(P))

  if (warn_missing && length(missing_genes) > 0) {
    warning(
      sprintf(
        "The following genes are not present in P and will be ignored: %s",
        paste(missing_genes, collapse = ", ")
      )
    )
  }

  if (length(present_genes) == 0) {
    stop("None of the requested genes are present in P.")
  }

  # Compute M_s = (1 / |G|) * sum_{g in G} P_gs
  axis_score <- colMeans(P[present_genes, , drop = FALSE], na.rm = na.rm)

  return(axis_score)
}

# ---------------------------
# Removal pathways (4 axes); gene sets are defined in R/common.R
# ---------------------------
rem_axis_list <- list(
  s1p = s1p_fin_axis,
  c1p = c1p_axis,
  sm  = sm_axis,
  gsl = gsl_axis
)

all_rem_targets <- unique(unlist(rem_axis_list))

viper_mat <- GetAssayData(metacell_obj, assay = "VIPER", layer = "data")

P_rem <- compute_target_pressure(
  regulon_obj = regulon_obj,
  viper_mat = viper_mat,
  targets = all_rem_targets
)

rem_pressure_scores <- lapply(rem_axis_list, function(genes) {
  compute_axis_score(P_rem, genes)
})

metacell_obj$Pressure_s1p <- rem_pressure_scores$s1p[colnames(metacell_obj)]
metacell_obj$Pressure_c1p <- rem_pressure_scores$c1p[colnames(metacell_obj)]
metacell_obj$Pressure_sm  <- rem_pressure_scores$sm[colnames(metacell_obj)]
metacell_obj$Pressure_gsl <- rem_pressure_scores$gsl[colnames(metacell_obj)]

out_rds <- file.path(out_dir, "metacell_obj_with_ceramide_axes_scores.rds")
saveRDS(metacell_obj, file = out_rds)
cat("Saved updated object to:", out_rds, "\n")
