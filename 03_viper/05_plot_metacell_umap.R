# =============================================================================
# 05. Figure 6a: metacell UMAP (SCT space) colored by condition.
#
# Input : <RESULTS_DIR>/03_viper/all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds
# Output: <RESULTS_DIR>/03_viper/figures/metacell_umap_condition_{square_no_legend,with_legend}.{svg,png}
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
})

source(file.path("R", "common.R"))

metacell_obj <- readRDS(file.path(DIR_VIPER, "all_metacell_hdWGCNA_with_VIPER_sct_umap_viper_umap.rds"))

# Set condition order
metacell_obj$condition <- factor(
  metacell_obj$condition,
  levels = CONDITION_LEVELS
)

# Base plot
p_base <- DimPlot(
  object = metacell_obj,
  reduction = "umap_sct_mc",
  group.by = "condition",
  label = FALSE,
  pt.size = 0.05
) +
  scale_color_manual(
    values = CONDITION_COLORS,
    breaks = CONDITION_LEVELS,
    drop = FALSE
  ) +
  labs(
    x = "UMAP_1",
    y = "UMAP_2",
    color = "Condition"
  ) +
  theme_classic(base_size = 7) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.35, "cm"),
    aspect.ratio = 1
  )

# 1) Square panel without legend
p_no_legend <- p_base +
  theme(
    legend.position = "none"
  )

# 2) With legend
p_with_legend <- p_base +
  theme(
    legend.position = "right"
  )

out_dir <- step_dir("03_viper", "figures")

save_plot_pair(p_no_legend, "metacell_umap_condition_square_no_legend", out_dir, w = 1.5, h = 1.5)
save_plot_pair(p_with_legend, "metacell_umap_condition_with_legend", out_dir, w = 2.5, h = 2.0)
