#!/usr/bin/env bash
# =============================================================================
# Run Cell Ranger ARC (cellranger-arc count) for the four samples in
# 00-1_cellranger_arc/samples.csv.
#
# Usage (from the repository root):
#   bash 00-1_cellranger_arc/run_cellranger_arc.sh            # all samples
#   bash 00-1_cellranger_arc/run_cellranger_arc.sh ct74       # one sample
#
# Settings (config/config.sh): CELLRANGER_ARC, ARC_REFERENCE, FASTQ_DIR,
# ARC_OUT_DIR, ARC_CORES, ARC_MEM_GB
# FASTQ layout: <FASTQ_DIR>/scRNA/<gex_fastq_dir>/*.fastq.gz
#               <FASTQ_DIR>/scATAC/<atac_fastq_dir>/*.fastq.gz
# Output:       <ARC_OUT_DIR>/<sample>/outs/  (use ARC_OUT_DIR as DATASET_DIR in config/config.R)
# =============================================================================
set -euo pipefail

source "${SNMULT_CONFIG_SH:-config/config.sh}"
SHEET="$(cd "$(dirname "$0")" && pwd)/samples.csv"
ONLY="${1:-}"

for v in CELLRANGER_ARC ARC_REFERENCE FASTQ_DIR ARC_OUT_DIR ARC_CORES ARC_MEM_GB; do
  [ -n "${!v:-}" ] || { echo "Set $v in config/config.sh" >&2; exit 2; }
done
[ -e "${ARC_REFERENCE}/reference.json" ] || { echo "Reference not found: ${ARC_REFERENCE}" >&2; exit 2; }

mkdir -p "${ARC_OUT_DIR}/libraries"
"${CELLRANGER_ARC}" --version

tail -n +2 "${SHEET}" | while IFS=, read -r sample condition description gex_dir gex_prefix atac_dir atac_prefix gex_index; do
  [ -z "${ONLY}" ] || [ "${ONLY}" = "${sample}" ] || continue

  if [ -e "${ARC_OUT_DIR}/${sample}" ]; then
    echo "[${sample}] ${ARC_OUT_DIR}/${sample} exists; skipping (remove it to rerun)"
    continue
  fi

  lib="${ARC_OUT_DIR}/libraries/${sample}.csv"
  printf 'fastqs,sample,library_type\n%s,%s,Gene Expression\n%s,%s,Chromatin Accessibility\n' \
    "${FASTQ_DIR}/scRNA/${gex_dir}" "${gex_prefix}" \
    "${FASTQ_DIR}/scATAC/${atac_dir}" "${atac_prefix}" > "${lib}"

  echo "[${sample}] ${condition}: cellranger-arc count"
  ( cd "${ARC_OUT_DIR}" && "${CELLRANGER_ARC}" count \
      --id="${sample}" \
      --reference="${ARC_REFERENCE}" \
      --libraries="${lib}" \
      --localcores="${ARC_CORES}" \
      --localmem="${ARC_MEM_GB}" < /dev/null )
done
