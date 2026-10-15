# =============================================================================
# 03. Metabolic diversity (Figure 6e, 6f).
#
#   (A) Per-metacell diversity: softmax (beta = 1) of the four axis pressures,
#       then Shannon entropy.
#   (B) Population diversity: 4D joint Shannon entropy per condition after
#       quantile binning of each axis (K = 6 bins).
#
# Input : <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores.rds
# Output: <RESULTS_DIR>/04_sphingolipid_axes/diversity_entropy_removal_axes/
# =============================================================================

source(file.path("R", "common.R"))

metacell_obj <- readRDS(file.path(DIR_AXES, "metacell_obj_with_ceramide_axes_scores.rds"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
})

# -----------------------------
# Output
# -----------------------------
sub_dir <- step_dir("04_sphingolipid_axes", "diversity_entropy_removal_axes")

# -----------------------------
# Data
# -----------------------------
md <- metacell_obj[[]]
need_cols <- c("condition", "Pressure_s1p", "Pressure_gsl", "Pressure_c1p", "Pressure_sm")
missing_cols <- setdiff(need_cols, colnames(md))
if (length(missing_cols) > 0) stop(sprintf("Missing meta.data columns: %s", paste(missing_cols, collapse = ", ")))

df <- data.frame(
  condition = factor(md$condition, levels = c("Control","Acute","Chronic","Return")),
  p_s1p = md$Pressure_s1p,
  p_gsl = md$Pressure_gsl,
  p_c1p = md$Pressure_c1p,
  p_sm  = md$Pressure_sm,
  stringsAsFactors = FALSE
)

ok <- is.finite(df$p_s1p) & is.finite(df$p_gsl) & is.finite(df$p_c1p) & is.finite(df$p_sm)
df <- df[ok, , drop = FALSE]

Pmat <- as.matrix(df[, c("p_s1p","p_gsl","p_c1p","p_sm"), drop = FALSE])
colnames(Pmat) <- c("S1P","GSL","C1P","SM")

# -----------------------------
# (A) Per-cell metabolic diversity (softmax -> entropy)
# -----------------------------
softmax_rows <- function(M, beta = 1.0) {
  X <- beta * M
  X <- X - apply(X, 1, max)
  EX <- exp(X)
  EX / rowSums(EX)
}

shannon_entropy_rows <- function(P, eps = 1e-12) {
  P2 <- pmax(P, eps)
  -rowSums(P2 * log(P2))
}

beta_pref <- 1.0
Pref <- softmax_rows(Pmat, beta = beta_pref)
df$metabolic_diversity <- shannon_entropy_rows(Pref)

condition_colors <- CONDITION_COLORS

p_cell <- ggplot(df, aes(x = condition, y = metabolic_diversity, fill = condition)) +
  geom_violin(scale = "width", trim = FALSE, linewidth = 0.2, color = "black") +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.2, fill = "white") +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  labs(
    title = sprintf("Per-cell metabolic diversity (softmax beta=%.2f)", beta_pref),
    x = NULL,
    y = "Metabolic diversity"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "none"
  )

ggsave(file.path(sub_dir, "diversity_cell_metabolic_violin.svg"),
       p_cell, width = 3.8, height = 2.6, units = "in", bg = "white")
ggsave(file.path(sub_dir, "diversity_cell_metabolic_violin.png"),
       p_cell, width = 3.8, height = 2.6, units = "in", dpi = 600, bg = "white")

# ---- minimal fix: safe summary table (same info, stable CSV)
summary_df <- data.frame(
  condition = levels(df$condition),
  mean = tapply(df$metabolic_diversity, df$condition, mean, na.rm = TRUE),
  median = tapply(df$metabolic_diversity, df$condition, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)
write.csv(
  summary_df,
  file.path(sub_dir, "diversity_cell_metabolic_summary.csv"),
  row.names = FALSE
)

# -----------------------------
# (B) Population diversity: 4D joint entropy by condition (K=6 quantile bins per axis)
# -----------------------------
K <- 6

make_quantile_breaks <- function(x, K) {
  qs <- quantile(x, probs = seq(0, 1, length.out = K + 1), na.rm = TRUE, names = FALSE)
  qs <- unique(qs)
  if (length(qs) < 3) stop("Too few unique quantile breaks; reduce K.")
  qs
}

joint_entropy_one <- function(ids) {
  tab <- table(ids)
  p <- as.numeric(tab) / sum(tab)
  H <- -sum(p * log(p))
  eff_states <- exp(H)
  c(H = H, eff_states = eff_states, n_bins_used = length(tab))
}

breaks_list <- lapply(as.data.frame(Pmat), make_quantile_breaks, K = K)

bin_mat <- sapply(seq_len(ncol(Pmat)), function(j) {
  cut(Pmat[, j], breaks = breaks_list[[j]], include.lowest = TRUE, labels = FALSE)
})
bin_id <- apply(bin_mat, 1, function(v) paste(v, collapse = "_"))

joint_list <- lapply(split(bin_id, df$condition), joint_entropy_one)
joint_out <- do.call(rbind, joint_list)

joint_df <- data.frame(
  condition = factor(rownames(joint_out), levels = c("Control","Acute","Chronic","Return")),
  population_diversity = joint_out[, "H"],
  eff_states = joint_out[, "eff_states"],
  n_bins_used = joint_out[, "n_bins_used"],
  K = K,
  stringsAsFactors = FALSE
)
rownames(joint_df) <- NULL

write.csv(joint_df, file.path(sub_dir, sprintf("diversity_population_joint4D_K%d.csv", K)), row.names = FALSE)

p_joint <- ggplot(joint_df, aes(x = condition, y = population_diversity, fill = condition)) +
  geom_col(color = "black", linewidth = 0.2) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_y_continuous(
    name = "Population diversity",
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = sprintf("Population diversity (4D joint entropy; quantile bins, K=%d)", K),
    x = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "none"
  )

ggsave(file.path(sub_dir, sprintf("diversity_population_joint4D_K%d_bar.svg", K)),
       p_joint, width = 3.8, height = 2.6, units = "in", bg = "white")
ggsave(file.path(sub_dir, sprintf("diversity_population_joint4D_K%d_bar.png", K)),
       p_joint, width = 3.8, height = 2.6, units = "in", dpi = 600, bg = "white")

cat("Saved metabolic (cell) and population (joint 4D) diversity plots to:", sub_dir, "\n")
