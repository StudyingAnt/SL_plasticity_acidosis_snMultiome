#!/usr/bin/env bash
# =============================================================================
# 00. Download the human transcription factor list used as ARACNe3 regulators
#     (Human TFs database v1.01; Lambert et al., Cell 2018) and check that it is
#     the same file used in the paper.
#
# Output: HUMAN_TF_FILE (config/config.sh), 1,639 gene symbols
#
# Usage (from the repository root):
#   bash 02_aracne3/00_download_human_tfs.sh
# =============================================================================
set -euo pipefail

source "${SNMULT_CONFIG_SH:-config/config.sh}"

URL="https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"
# SHA-256 of the file used in the paper
EXPECTED_SHA256="52fc1dab0d02b3a28b3bbe45232fceb2f85ebafd23631f621c5efc8bc91767c2"

mkdir -p "$(dirname "${HUMAN_TF_FILE}")"

if command -v wget >/dev/null 2>&1; then
  wget -q -O "${HUMAN_TF_FILE}" "${URL}"
else
  curl -fsSL -o "${HUMAN_TF_FILE}" "${URL}"
fi

actual=$(sha256sum "${HUMAN_TF_FILE}" | awk '{print $1}')
echo "Saved: ${HUMAN_TF_FILE} ($(grep -c . "${HUMAN_TF_FILE}") TFs)"

if [ "${actual}" = "${EXPECTED_SHA256}" ]; then
  echo "Checksum OK: identical to the file used in the paper."
else
  echo "WARNING: checksum differs from the file used in the paper." >&2
  echo "  expected ${EXPECTED_SHA256}" >&2
  echo "  actual   ${actual}" >&2
  exit 1
fi
