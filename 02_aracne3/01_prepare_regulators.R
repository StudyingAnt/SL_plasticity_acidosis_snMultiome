# =============================================================================
# 01. Prepare the ARACNe3 regulator list: human TFs (Human TFs database v1.01)
#     that are present in the metacell CPM matrix.
#
# Input : HUMAN_TF_FILE (config/config.R)
#         <RESULTS_DIR>/01_metacells/metacell_expression_CPM.tsv
# Output: <RESULTS_DIR>/02_aracne3/regulators_in_matrix.txt
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(data.table)
})

source(file.path("R", "common.R"))

out_dir <- step_dir("02_aracne3")

cpm_mat <- read_cpm_matrix()

# Human TF file(TF gene names) was manually downloaded from: https://humantfs.ccbr.utoronto.ca/download.php
human_tf_file <- HUMAN_TF_FILE

tf_symbols <- readLines(human_tf_file, warn = FALSE)
tf_symbols <- unique(trimws(tf_symbols))
tf_symbols <- tf_symbols[nzchar(tf_symbols)]

regs_keep <- intersect(tf_symbols, rownames(cpm_mat))

regulator_file_out <- file.path(out_dir, "regulators_in_matrix.txt")
writeLines(regs_keep, con = regulator_file_out)

cat(glue("[TF] Total TFs in source file: {length(tf_symbols)}\n"))
cat(glue("[TF] TFs present in cpm_mat: {length(regs_keep)}\n"))
cat(glue("[TF] Saved regulator list: {regulator_file_out}\n"))
