# =============================================================================
# Shared setup for 00-2_seurat_preprocessing (sourced by 01-07).
# Loads packages, sample table, and hg38 gene annotation.
# =============================================================================

suppressMessages({
  library(Seurat)
  library(Signac)

  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)

  library(dplyr)
  library(tibble)
  library(fs)
  library(glue)

  library(jsonlite)
})

source(file.path("R", "common.R"))

# I/O
seurat_dir <- SEURAT_DIR
qc_plots_dir <- step_dir("00-2_seurat_preprocessing", "qc_plots")
macs_out_root <- step_dir("00-2_seurat_preprocessing", "macs3_out")
dir_create(seurat_dir)

cat("DATASET_DIR:", DATASET_DIR, "\n")
cat("seurat_dir:", seurat_dir, "\n")

# ----------------------------
# Samples / mapping
# ----------------------------
samples <- as_tibble(SAMPLES) |>
  mutate(
    ess_dir = file.path(DATASET_DIR, sample),
    cnts_h5 = file.path(ess_dir, "filtered_feature_bc_matrix.h5"),
    frags   = file.path(ess_dir, "atac_fragments.tsv.gz")
  )

stopifnot(all(file_exists(samples$cnts_h5)))
stopifnot(all(file_exists(samples$frags)))

# ----------------------------
# Annotation (hg38) from EnsDb
# ----------------------------
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
