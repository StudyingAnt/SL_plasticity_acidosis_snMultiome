# =============================================================================
# 06. ACDC-style grid search clustering (Leiden; k = 11..101, resolution =
#     0.1..1.9) on pca.sct and peaks.lsi, selecting the parameters with the
#     highest mean per-cluster silhouette width.
# Usage: Rscript 00-2_seurat_preprocessing/06_clustering_acdc.R \
#          -i <SEURAT_DIR>/all_qc_noMT_macs3_trf.rds
# Output: <SEURAT_DIR>/<name>_qc_noMT_macs3_trf_clust.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))
suppressMessages({
  library(cluster)
  library(future)
})

# Same limit as in 05_transformations.R; FindNeighbors exports ~700 MiB to
# future workers for the merged object (size limit only, no effect on results)
options(future.globals.maxSize = 10 * 1024^3)

# Parse CLI args (base R)
args <- commandArgs(trailingOnly = TRUE)

input_rds <- NA_character_
if (length(args) > 0) {
  i_pos <- which(args %in% c("-i", "--input"))
  if (length(i_pos) > 0) {
    if (i_pos[1] == length(args)) stop("'-i/--input' requires a path.")
    input_rds <- args[i_pos[1] + 1]
  }
}

if (is.na(input_rds)) {
  stop("Missing input. Use: Rscript 00-2_seurat_preprocessing/06_clustering_acdc.R -i /path/to/sample.rds")
}

obj_file <- input_rds
if (!file.exists(obj_file)) {
  stop(glue("Input RDS not found: {obj_file}"))
}

# ----------------------------
# Load object
# ----------------------------
cat(glue("Loading SeuratObject from: {obj_file}\n"))
obj <- readRDS(obj_file)

# Infer obj_name: split filename by '_' and take the first token
bn <- basename(obj_file)
obj_name <- strsplit(bn, "_", fixed = TRUE)[[1]][1]

cat(glue(
  "[{obj_name}] SeuratObject loaded.\n",
  "  Assays: {paste(Seurat::Assays(obj), collapse = ', ')}\n",
  "  Cells: {ncol(obj)}\n\n"
))

# SCT
# --------
# Settings
# --------
KNN_GRID <- seq(11, 101, by = 10)
RES_GRID <- seq(0.1, 1.9, by = 0.2)
DIMS_SCT <- 1:50

# --------
# Checks
# --------
if (!("pca.sct" %in% Reductions(obj))) {
  stop(glue("[{obj_name}] Missing reduction 'pca.sct'."))
}

cat(glue("\n[{obj_name}] [SCT] Precomputing distance matrix (Euclidean on pca.sct dims {min(DIMS_SCT)}:{max(DIMS_SCT)})...\n"))
emb_sct <- Embeddings(obj, "pca.sct")[, DIMS_SCT, drop = FALSE]

# dist object (stores n*(n-1)/2 distances)
dist_sct <- dist(emb_sct)

cat(glue("[{obj_name}] [SCT] Distance computed. Starting grid search...\n"))

# Result table (like GS.history)
sct_hist <- tibble::tibble(
  knn = integer(0),
  resolution = numeric(0),
  sil_avg = numeric(0), # ACDC-like: mean of per-cluster mean silhouettes (unweighted)
  sil_overall = numeric(0), # overall mean silhouette (cell-weighted)
  n_clust = integer(0)
)

# Use fixed tmp graph/cluster names to avoid object bloat
tmp_nn <- "sct_tmp.nn"
tmp_snn <- "sct_tmp.snn"
tmp_col <- "sct_tmp_clusters"

for (K in KNN_GRID) {
  cat(glue("[{obj_name}] [SCT] KNN={K}: FindNeighbors -> graphs ({tmp_nn}, {tmp_snn})\n"))

  obj <- FindNeighbors(
    object = obj,
    reduction = "pca.sct",
    dims = DIMS_SCT,
    k.param = K,
    graph.name = c(tmp_nn, tmp_snn),
    verbose = FALSE
  )

  for (R in RES_GRID) {
    obj <- FindClusters(
      object = obj,
      graph.name = tmp_snn,
      resolution = R,
      algorithm = 4, # Leiden
      cluster.name = tmp_col, # overwrite same column each time
      verbose = FALSE
    )

    cl <- obj[[tmp_col, drop = TRUE]]
    ncl <- length(unique(cl))
    if (ncl <= 1) next

    sil <- cluster::silhouette(as.integer(as.factor(cl)), dist_sct)

    sil_overall <- mean(sil[, "sil_width"])
    sil_by_cl <- tapply(sil[, "sil_width"], cl, mean)
    sil_avg <- mean(sil_by_cl)

    sct_hist <- dplyr::bind_rows(
      sct_hist,
      tibble::tibble(knn = K, resolution = R, sil_avg = sil_avg, sil_overall = sil_overall, n_clust = ncl)
    )
  }
}

if (nrow(sct_hist) == 0) {
  stop(glue("[{obj_name}] [SCT] No valid clustering (all n_clust <= 1)."))
}

# Pick best (tie -> first)
i_best <- which(sct_hist$sil_avg == max(sct_hist$sil_avg, na.rm = TRUE))[1]
best <- sct_hist[i_best, , drop = FALSE]

cat(glue(
  "[{obj_name}] [SCT] Best found:\n",
  "  knn={best$knn}, res={best$resolution}\n",
  "  sil_avg={signif(best$sil_avg, 4)}, sil_overall={signif(best$sil_overall, 4)}\n",
  "  n_clust={best$n_clust}\n"
))

# Store history + best params
obj[["SCT"]]@misc$GS.history <- sct_hist
obj[["SCT"]]@misc$GS.opt_params <- best[, c("knn", "resolution"), drop = FALSE]

# Final run (keep separate final graph names if you want reproducibility)
obj <- FindNeighbors(
  object = obj,
  reduction = "pca.sct",
  dims = DIMS_SCT,
  k.param = best$knn,
  graph.name = c("sct_opt.nn", "sct_opt.snn"),
  verbose = FALSE
)
obj <- FindClusters(
  object = obj,
  graph.name = "sct_opt.snn",
  resolution = best$resolution,
  algorithm = 4,
  cluster.name = "acdc_sct",
  verbose = FALSE
)

# Cleanup tmp (optional but recommended)
obj@graphs[[tmp_nn]] <- NULL
obj@graphs[[tmp_snn]] <- NULL
obj@meta.data[[tmp_col]] <- NULL

cat(glue("[{obj_name}] [SCT] Final clusters saved in meta.data$acdc_sct\n"))

# peaks
KNN_GRID <- seq(11, 101, by = 10)
RES_GRID <- seq(0.1, 1.9, by = 0.2)
DIMS_LSI <- 2:51 # skip LSI_1

if (!("peaks.lsi" %in% Reductions(obj))) {
  stop(glue("[{obj_name}] Missing reduction 'peaks.lsi'."))
}

cat(glue("\n[{obj_name}] [peaks] Precomputing distance matrix (Euclidean on peaks.lsi dims {min(DIMS_LSI)}:{max(DIMS_LSI)})...\n"))
emb_peaks <- Embeddings(obj, "peaks.lsi")[, DIMS_LSI, drop = FALSE]
dist_peaks <- dist(emb_peaks)

cat(glue("[{obj_name}] [peaks] Distance computed. Starting grid search...\n"))

peaks_hist <- tibble::tibble(
  knn = integer(0),
  resolution = numeric(0),
  sil_avg = numeric(0),
  sil_overall = numeric(0),
  n_clust = integer(0)
)

tmp_nn <- "peaks_tmp.nn"
tmp_snn <- "peaks_tmp.snn"
tmp_col <- "peaks_tmp_clusters"

for (K in KNN_GRID) {
  cat(glue("[{obj_name}] [peaks] KNN={K}: FindNeighbors -> graphs ({tmp_nn}, {tmp_snn})\n"))

  obj <- FindNeighbors(
    object = obj,
    reduction = "peaks.lsi",
    dims = DIMS_LSI,
    k.param = K,
    graph.name = c(tmp_nn, tmp_snn),
    verbose = FALSE
  )

  for (R in RES_GRID) {
    obj <- FindClusters(
      object = obj,
      graph.name = tmp_snn,
      resolution = R,
      algorithm = 4,
      cluster.name = tmp_col,
      verbose = FALSE
    )

    cl <- obj[[tmp_col, drop = TRUE]]
    ncl <- length(unique(cl))
    if (ncl <= 1) next

    sil <- cluster::silhouette(as.integer(as.factor(cl)), dist_peaks)

    sil_overall <- mean(sil[, "sil_width"])
    sil_by_cl <- tapply(sil[, "sil_width"], cl, mean)
    sil_avg <- mean(sil_by_cl)

    peaks_hist <- dplyr::bind_rows(
      peaks_hist,
      tibble::tibble(knn = K, resolution = R, sil_avg = sil_avg, sil_overall = sil_overall, n_clust = ncl)
    )
  }
}

if (nrow(peaks_hist) == 0) {
  stop(glue("[{obj_name}] [peaks] No valid clustering (all n_clust <= 1)."))
}

i_best <- which(peaks_hist$sil_avg == max(peaks_hist$sil_avg, na.rm = TRUE))[1]
best <- peaks_hist[i_best, , drop = FALSE]

cat(glue(
  "[{obj_name}] [peaks] Best found:\n",
  "  knn={best$knn}, res={best$resolution}\n",
  "  sil_avg={signif(best$sil_avg, 4)}, sil_overall={signif(best$sil_overall, 4)}\n",
  "  n_clust={best$n_clust}\n"
))

obj[["peaks"]]@misc$GS.history <- peaks_hist
obj[["peaks"]]@misc$GS.opt_params <- best[, c("knn", "resolution"), drop = FALSE]

obj <- FindNeighbors(
  object = obj,
  reduction = "peaks.lsi",
  dims = DIMS_LSI,
  k.param = best$knn,
  graph.name = c("peaks_opt.nn", "peaks_opt.snn"),
  verbose = FALSE
)
obj <- FindClusters(
  object = obj,
  graph.name = "peaks_opt.snn",
  resolution = best$resolution,
  algorithm = 4,
  cluster.name = "acdc_peaks",
  verbose = FALSE
)

# Cleanup tmp
obj@graphs[[tmp_nn]] <- NULL
obj@graphs[[tmp_snn]] <- NULL
obj@meta.data[[tmp_col]] <- NULL

cat(glue("[{obj_name}] [peaks] Final clusters saved in meta.data$acdc_peaks\n"))

# ----------------------------
# Save output
# ----------------------------
out_rds <- file.path(seurat_dir, paste0(obj_name, "_qc_noMT_macs3_trf_clust.rds"))
cat(glue("[{obj_name}] Saving clustered SeuratObject...\n"))
saveRDS(obj, out_rds)
cat(glue("[{obj_name}] Saved in: {out_rds}\n\n"))

cat(glue("[DONE] ACDC GridSearch clustering + saving completed for all samples and merged object.\n"))
