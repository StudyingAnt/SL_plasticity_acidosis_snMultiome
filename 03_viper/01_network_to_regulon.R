# =============================================================================
# 01. Convert the ARACNe3 consolidated network to 3-column format and build a
#     VIPER regulon object with aracne2regulon.
#
# Input : <RESULTS_DIR>/02_aracne3/consolidated-net_mcf7_metacell.tsv
#         <RESULTS_DIR>/01_metacells/metacell_expression_CPM.tsv
# Output: <RESULTS_DIR>/03_viper/consolidated-net_mcf7_metacell_for_viper_3col.tsv
#         <RESULTS_DIR>/03_viper/mcf7_metacell_regulon_viper.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(viper)
  library(glue)
})

source(file.path("R", "common.R"))

aracne_out_dir <- step_dir("03_viper")

# -----------------------------------------------------------------------------
# Change ARACNe3 network file for VIPER
# -----------------------------------------------------------------------------
aracne3_net_file <- file.path(DIR_ARACNE3, "consolidated-net_mcf7_metacell.tsv")
viper_net_file <- file.path(aracne_out_dir, "consolidated-net_mcf7_metacell_for_viper_3col.tsv")

net_dt <- fread(aracne3_net_file)

cat(glue("[VIPER] Loaded ARACNe3 network: {nrow(net_dt)} edges, {ncol(net_dt)} columns\n"))
cat(glue("[VIPER] Columns: {paste(colnames(net_dt), collapse = ', ')}\n"))

if (ncol(net_dt) < 3) {
  stop("ARACNe3 network file has fewer than 3 columns.")
}

net3 <- net_dt[, 1:3]

# Remove possible NA rows
net3 <- net3[!is.na(net3[[1]]) & !is.na(net3[[2]]) & !is.na(net3[[3]])]

# Write without header for maximum compatibility
fwrite(
  net3,
  file = viper_net_file,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)

cat(glue("[VIPER] Saved 3-column network file: {viper_net_file}\n"))

# -----------------------------------------------------------------------------
# Load expression matrix
# -----------------------------------------------------------------------------
expr_mat <- read_cpm_matrix()

cat(glue("[VIPER] Expression matrix: {nrow(expr_mat)} genes x {ncol(expr_mat)} metacells\n"))

# -----------------------------------------------------------------------------
# Make regulon object
# -----------------------------------------------------------------------------
regulon_obj <- aracne2regulon(
  afile = viper_net_file,
  eset = expr_mat,
  format = "3col",
  verbose = TRUE
)

regulon_rds <- file.path(aracne_out_dir, "mcf7_metacell_regulon_viper.rds")
saveRDS(regulon_obj, regulon_rds)

cat(glue("[VIPER] Saved regulon object: {regulon_rds}\n"))
print(regulon_obj)
