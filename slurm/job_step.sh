#!/usr/bin/env bash
# =============================================================================
# Body of a SLURM job that runs one pipeline step. Submitted by slurm/submit.sh;
# not meant to be called directly.
#   $1     step script (path relative to the repository root)
#   $2...  arguments passed to the step script
# =============================================================================
set -eo pipefail

# slurm/submit.sh submits from the repository root
cd "${SLURM_SUBMIT_DIR}"

source config/config.sh

# R environment (conda R + project R library)
if [ -n "${CONDA_MODULE:-}" ]; then
  module load "${CONDA_MODULE}"
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
export R_LIBS="${R_LIBRARY}"

step="$1"
shift

echo "=== job ${SLURM_JOB_ID} on $(hostname) | $(date)"
echo "=== repo commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "=== step: ${step} $*"
echo "=== CPUs: ${SLURM_CPUS_PER_TASK:-?}  memory: ${SLURM_MEM_PER_NODE:-?} MB"

case "${step}" in
  *.R)  Rscript "${step}" "$@" ;;
  *.sh) bash "${step}" "$@" ;;
  *)    echo "Unknown step type: ${step}" >&2; exit 1 ;;
esac

echo "=== finished $(date)"
if [[ "${step}" == *.R ]]; then
  Rscript -e '
pk <- c("Seurat", "SeuratObject", "Signac", "sctransform", "glmGamPoi", "chromVAR", "motifmatchr",
        "JASPAR2024", "TFBSTools", "hdWGCNA", "viper", "genefilter", "data.table", "ggplot2", "Matrix")
cat(R.version.string, "\n")
for (p in pk) cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")))
'
fi
