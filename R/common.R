# =============================================================================
# Shared settings sourced by every R script in this repository.
# Run all scripts from the repository root.
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
})

# SNMULT_CONFIG_R can point to another config file (e.g. for a test run)
cfg_file <- Sys.getenv("SNMULT_CONFIG_R", file.path("config", "config.R"))
if (!file.exists(cfg_file)) {
  stop(
    cfg_file, " not found. Copy config/config.example.R to config/config.R, ",
    "edit the paths, and run scripts from the repository root."
  )
}
cat("Config:", cfg_file, "\n")
source(cfg_file)

# ----------------------------
# Samples / conditions
# ----------------------------
SAMPLES <- data.frame(
  sample = c("ct74", "ct65", "aa65", "aa74"),
  condition = c("Control", "Acute", "Chronic", "Return"),
  stringsAsFactors = FALSE
)

CONDITION_LEVELS <- c("Control", "Acute", "Chronic", "Return")

CONDITION_COLORS <- c(
  "Control" = "#2166AC",
  "Acute"   = "#C51B7D",
  "Chronic" = "#B30000",
  "Return"  = "#FC8D59"
)

# ----------------------------
# Ceramide clearance axes (KEGG hsa00600, manually curated)
# ----------------------------
s1p_fin_axis <- c("ACER1", "ACER2", "ASAH1", "ASAH2", "SPHK1", "SPHK2")
c1p_axis <- c("CERK")
sm_axis <- c("SGMS1", "SGMS2")
gsl_axis <- c("UGCG", "UGT8", "B4GALT4", "B4GALT5")

# ----------------------------
# Output directories (one per analysis step)
# ----------------------------
step_dir <- function(...) {
  d <- file.path(RESULTS_DIR, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

DIR_METACELLS <- file.path(RESULTS_DIR, "01_metacells")
DIR_ARACNE3   <- file.path(RESULTS_DIR, "02_aracne3")
DIR_VIPER     <- file.path(RESULTS_DIR, "03_viper")
DIR_AXES      <- file.path(RESULTS_DIR, "04_sphingolipid_axes")
DIR_REGS      <- file.path(RESULTS_DIR, "05_axis_regulators")

# ----------------------------
# Plot helper: save the same ggplot as SVG and PNG
# ----------------------------
save_plot_pair <- function(p, out_prefix, out_dir, w, h, dpi = 600) {
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, ".svg")),
    plot = p, width = w, height = h, units = "in", bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(out_prefix, ".png")),
    plot = p, width = w, height = h, units = "in", dpi = dpi, bg = "white"
  )
}

# ----------------------------
# Read the metacell CPM matrix (genes x metacells) written by 01_metacells
# ----------------------------
read_cpm_matrix <- function(file = file.path(DIR_METACELLS, "metacell_expression_CPM.tsv")) {
  expr_dt <- data.table::fread(file)
  expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
  rownames(expr_mat) <- expr_dt[[1]]
  mode(expr_mat) <- "numeric"
  expr_mat
}
