#!/usr/bin/env bash
# =============================================================================
# Local configuration for shell steps and SLURM (template)
#
# Copy this file to config/config.sh and edit. config/config.sh is ignored by git.
# RESULTS_DIR must point to the same directory as RESULTS_DIR in config/config.R.
# =============================================================================

RESULTS_DIR="/path/to/work/results"

# Human TF list (same path as HUMAN_TF_FILE in config/config.R)
HUMAN_TF_FILE="/path/to/TF_names_v_1.01.txt"

# ARACNe3 executable (https://github.com/califano-lab/ARACNe3)
ARACNE3_BIN="ARACNe3_app_release"

ARACNE3_THREADS=8

# -----------------------------------------------------------------------------
# SLURM (used by slurm/submit.sh and slurm/job_step.sh)
# -----------------------------------------------------------------------------
# Environment module that provides conda (leave empty if conda is already on PATH)
CONDA_MODULE=""

# conda environment with R (path or name)
CONDA_ENV="/path/to/envs/r452"

# R package library used for the analysis (e.g. a restored renv library)
R_LIBRARY="/path/to/r/library"

# SLURM partition and log directory (keep logs outside the repository)
SLURM_PARTITION="partition-name"
LOG_DIR="/path/to/work/logs"

# Extra sbatch options, e.g. "--account=xxx --mail-user=me@example.org --mail-type=END,FAIL"
SLURM_EXTRA_ARGS=""
