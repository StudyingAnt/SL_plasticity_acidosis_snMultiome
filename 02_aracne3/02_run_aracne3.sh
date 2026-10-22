#!/usr/bin/env bash
# =============================================================================
# 02. Infer the TF regulatory network from metacell CPM with ARACNe3.
#
# Input : <RESULTS_DIR>/01_metacells/metacell_expression_CPM.tsv
#         <RESULTS_DIR>/02_aracne3/regulators_in_matrix.txt
# Output: <RESULTS_DIR>/02_aracne3/consolidated-net_mcf7_metacell.tsv (+ subnets, logs)
#
# Usage (from the repository root):
#   bash 02_aracne3/02_run_aracne3.sh
# =============================================================================
set -euo pipefail

source "${SNMULT_CONFIG_SH:-config/config.sh}"

ARACNe3_OUT="${RESULTS_DIR}/02_aracne3"
mkdir -p "${ARACNe3_OUT}"

# ARACNe3 caches its null model in ./.ARACNe3_cached/ of the working directory.
# A run that loads this cache skips the random draws used to build the null
# model, so the bootstrap subsamples (and the network) differ from a fresh run
# with the same --seed. Run inside the output directory and start without a
# cache so every run matches the run used in the paper.
cd "${ARACNe3_OUT}"
rm -rf .ARACNe3_cached

"${ARACNE3_BIN}" \
  -e "${RESULTS_DIR}/01_metacells/metacell_expression_CPM.tsv" \
  -r "${ARACNe3_OUT}/regulators_in_matrix.txt" \
  -o "${ARACNe3_OUT}" \
  -x 30 \
  --alpha 0.05 \
  --threads "${ARACNE3_THREADS}" \
  --seed 1234 \
  --runid mcf7_metacell
