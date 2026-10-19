#!/usr/bin/env bash
# =============================================================================
# Submit one pipeline step to SLURM.
#
# Usage (from anywhere):
#   bash slurm/submit.sh <step script> [script arguments]
# Examples:
#   bash slurm/submit.sh 00-2_seurat_preprocessing/01_make_objects.R
#   bash slurm/submit.sh 00-2_seurat_preprocessing/06_clustering_acdc.R -i /path/all_qc_noMT_macs3_trf.rds
#
# Resources below are defaults per step; override with environment variables,
# e.g.  MEM=200G TIME=12:00:00 bash slurm/submit.sh <step>
# Partition, log directory and extra sbatch options come from config/config.sh.
#
# To use other config files (e.g. a test run with its own results directory):
#   SNMULT_CONFIG_R=config/config.test.R SNMULT_CONFIG_SH=config/config.test.sh \
#     JOB_TAG=test bash slurm/submit.sh <step>
# =============================================================================
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO}"

CONFIG_SH="${SNMULT_CONFIG_SH:-config/config.sh}"
if [ ! -f "${CONFIG_SH}" ]; then
  echo "${CONFIG_SH} not found. Copy config/config.example.sh to config/config.sh and edit it." >&2
  exit 1
fi
source "${CONFIG_SH}"

if [ $# -lt 1 ]; then
  echo "Usage: bash slurm/submit.sh <step script> [arguments]" >&2
  exit 1
fi
step="$1"
if [ ! -f "${step}" ]; then
  echo "Step script not found: ${step}" >&2
  exit 1
fi

# cpus  mem  time
case "${step}" in
  00-2_seurat_preprocessing/01_make_objects.R)       res="4  64G  08:00:00" ;;
  00-2_seurat_preprocessing/02_qc_objects.R)         res="2  64G  04:00:00" ;;
  00-2_seurat_preprocessing/03_call_peaks.R)         res="8  128G 24:00:00" ;;
  00-2_seurat_preprocessing/04_merge_objects.R)      res="2  128G 04:00:00" ;;
  00-2_seurat_preprocessing/05_transformations.R)    res="8  240G 24:00:00" ;;
  00-2_seurat_preprocessing/06_clustering_acdc.R)    res="4  200G 24:00:00" ;;
  00-2_seurat_preprocessing/07_run_chromvar.R)       res="8  200G 24:00:00" ;;
  01_metacells/01_build_metacells.R)                 res="4  200G 12:00:00" ;;
  02_aracne3/01_prepare_regulators.R)                res="1  16G  01:00:00" ;;
  02_aracne3/02_run_aracne3.sh)                      res="8  64G  24:00:00" ;;
  03_viper/04_attach_chromvar.R)                     res="2  200G 08:00:00" ;;
  *)                                                 res="4  64G  08:00:00" ;;
esac
read -r d_cpus d_mem d_time <<< "${res}"
CPUS="${CPUS:-${d_cpus}}"
MEM="${MEM:-${d_mem}}"
TIME="${TIME:-${d_time}}"

name="snmult${JOB_TAG:+_${JOB_TAG}}_$(basename "${step}" | sed 's/\.[^.]*$//')"
mkdir -p "${LOG_DIR}"

# shellcheck disable=SC2086
sbatch \
  --job-name="${name}" \
  --partition="${SLURM_PARTITION}" \
  --nodes=1 --ntasks=1 \
  --cpus-per-task="${CPUS}" \
  --mem="${MEM}" \
  --time="${TIME}" \
  --output="${LOG_DIR}/%x.%j.out" \
  --error="${LOG_DIR}/%x.%j.err" \
  ${SLURM_EXTRA_ARGS:-} \
  slurm/job_step.sh "$@"
