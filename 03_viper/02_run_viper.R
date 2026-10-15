# =============================================================================
# 02. Run VIPER (rank method, minsize 25) on the metacell CPM matrix and attach
#     the activity matrix to the metacell object as assay "VIPER".
#
# Input : <RESULTS_DIR>/01_metacells/all_metacell_hdWGCNA.rds
#         <RESULTS_DIR>/01_metacells/metacell_expression_CPM.tsv
#         <RESULTS_DIR>/03_viper/mcf7_metacell_regulon_viper.rds
# Output: <RESULTS_DIR>/03_viper/mcf7_metacell_viper_rank.{rds,tsv}
#         <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER.rds
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(Seurat)
  library(SeuratObject)
  library(viper)
  library(data.table)
})

source(file.path("R", "common.R"))

aracne_out_dir <- step_dir("03_viper")

# Load saved objects
metacell_obj <- readRDS(file.path(DIR_METACELLS, "all_metacell_hdWGCNA.rds"))
regulon_obj <- readRDS(file.path(aracne_out_dir, "mcf7_metacell_regulon_viper.rds"))

# Load expression matrix used for VIPER
expr_mat <- read_cpm_matrix()

# Run VIPER
viper_mat <- viper(
  eset = expr_mat,
  regulon = regulon_obj,
  method = "rank",
  minsize = 25,
  nes = TRUE,
  pleiotropy = FALSE,
  cores = 1,
  verbose = TRUE
)

# Save raw VIPER matrix immediately
saveRDS(
  viper_mat,
  file.path(aracne_out_dir, "mcf7_metacell_viper_rank.rds")
)
fwrite(
  data.table(regulator = rownames(viper_mat), as.data.frame(viper_mat)),
  file.path(aracne_out_dir, "mcf7_metacell_viper_rank.tsv"),
  sep = "\t",
  quote = FALSE
)

cat(glue("[VIPER] Saved activity matrix.\n"))
cat(glue("[VIPER] Matrix size: {nrow(viper_mat)} regulators x {ncol(viper_mat)} metacells\n"))

# Attach VIPER assay to Seurat object
metacell_obj[["VIPER"]] <- CreateAssayObject(data = viper_mat)

# Save object with VIPER assay
saveRDS(
  metacell_obj,
  file.path(aracne_out_dir, "all_metacell_hdWGCNA_with_VIPER.rds")
)

cat(glue("[VIPER] Saved Seurat object with VIPER assay.\n"))
cat(glue("Assays: {paste(names(metacell_obj@assays), collapse = ', ')}\n"))
