# =============================================================================
# Local configuration (template)
#
# Copy this file to config/config.R and edit the paths for your machine.
# config/config.R is ignored by git.
# All R scripts are meant to be run from the repository root, e.g.
#   Rscript 01_metacells/01_build_metacells.R
# =============================================================================

# Cell Ranger ARC outputs (see 00-1_cellranger_arc; same path as ARC_OUT_DIR).
# Expected layout: <DATASET_DIR>/<sample>/outs/filtered_feature_bc_matrix.h5
#                  <DATASET_DIR>/<sample>/outs/atac_fragments.tsv.gz
# (files placed directly in <DATASET_DIR>/<sample>/ also work)
DATASET_DIR <- "/path/to/cellranger_arc_outputs"

# Intermediate Seurat objects written by 00-2_seurat_preprocessing
SEURAT_DIR <- "/path/to/work/seurat_objects"

# All downstream results (tables, RDS objects, figures)
RESULTS_DIR <- "/path/to/work/results"

# MACS3 executable used by Signac::CallPeaks (00-2_seurat_preprocessing/03_call_peaks.R)
MACS3_PATH <- "/path/to/macs3"

# Human TF list (Human TFs database v1.01, TF_names_v_1.01.txt), see 02_aracne3/README.md
HUMAN_TF_FILE <- "/path/to/TF_names_v_1.01.txt"
