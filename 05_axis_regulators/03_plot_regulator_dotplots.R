# =============================================================================
# 03. Figure 6f-i: candidate regulators per axis.
#     msVIPER padj < 0.05 (full-transcriptome universe); top 20 regulators per
#     axis by mean -log10(padj) across conditions.
#     color = NES, size = axis gene linkage, alpha = -log10(motif padj).
#
# Input : <RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/final_candidates_topbottom20.tsv
# Output: <RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/viz_allgenes_regulon_motifalpha_padj005/
#           dotplot_padj005_nes_linkage_motifalpha_<axis>.{svg,png}
#           plot_dt_allgenes_regulon_motifalpha_padj005.tsv
# =============================================================================

source(file.path("R", "common.R"))

msviper_tb_dir <- file.path(DIR_REGS, "msviper_topbottom20_by_axis")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

plot_dir <- file.path(msviper_tb_dir, "viz_allgenes_regulon_motifalpha_padj005")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

final_file <- file.path(msviper_tb_dir, "final_candidates_topbottom20.tsv")
if (!file.exists(final_file)) stop("final_candidates_topbottom20.tsv not found.")

plot_dt <- fread(final_file)
plot_dt <- plot_dt[universe == "all_genes"]

plot_dt[, condition := tolower(condition)]
plot_dt[, axis := toupper(axis)]

plot_dt_filt <- plot_dt[padj < 0.05]

plot_dt_filt[, condition_label := factor(
  condition,
  levels = c("control", "acute", "chronic", "return"),
  labels = c("Control", "Acute", "Chronic", "Return")
)]

plot_dt_filt[is.na(motif_best_padj), motif_best_padj := 1]
plot_dt_filt[, motif_logp := -log10(motif_best_padj + 1e-300)]

topM <- 20

nes_lim <- max(abs(plot_dt_filt$nes), na.rm = TRUE)
linkage_lim <- range(plot_dt_filt$sum_abs_w, na.rm = TRUE, finite = TRUE)
motif_alpha_lim <- range(plot_dt_filt$motif_logp, na.rm = TRUE, finite = TRUE)

for (ax in unique(plot_dt_filt$axis)) {
  regs <- plot_dt_filt[
    axis == ax,
    .(mean_msviper_logp = mean(-log10(padj + 1e-300), na.rm = TRUE)),
    by = regulator
  ][order(-mean_msviper_logp)][1:min(.N, topM), regulator]

  d <- plot_dt_filt[axis == ax & regulator %in% regs]
  if (nrow(d) == 0) next

  reg_order <- d[, .(mean_msviper_logp = mean(-log10(padj + 1e-300), na.rm = TRUE)), by = regulator][
    order(mean_msviper_logp)
  ][, regulator]

  d[, regulator := factor(regulator, levels = reg_order)]

  p <- ggplot(d, aes(x = condition_label, y = regulator)) +
    geom_point(aes(
      color = nes,
      size = sum_abs_w,
      alpha = motif_logp
    )) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-nes_lim, nes_lim),
      name = "NES"
    ) +
    scale_size_continuous(
      limits = linkage_lim,
      name = "Axis gene linkage"
    ) +
    scale_alpha_continuous(
      limits = motif_alpha_lim,
      range = c(0.25, 1),
      name = "-log10(motif padj)"
    ) +
    guides(
      color = guide_colorbar(
        order = 1,
        barheight = unit(18, "mm"),
        barwidth = unit(4, "mm")
      ),
      size = guide_legend(
        order = 2,
        keyheight = unit(3.5, "mm"),
        keywidth = unit(3.5, "mm")
      ),
      alpha = guide_legend(
        order = 3,
        keyheight = unit(3.5, "mm"),
        keywidth = unit(3.5, "mm")
      )
    ) +
    labs(
      title = paste0(ax, " Axis"),
      x = NULL,
      y = NULL
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 6),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(3.5, "mm"),
      legend.spacing.y = unit(1.5, "mm")
    )

  save_plot_pair(
    p,
    out_prefix = paste0("dotplot_padj005_nes_linkage_motifalpha_", ax),
    out_dir = plot_dir,
    w = 2.7,
    h = 3.3,
    dpi = 600
  )
}

fwrite(
  plot_dt_filt,
  file.path(plot_dir, "plot_dt_allgenes_regulon_motifalpha_padj005.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("Saved filtered dot plots to:", plot_dir, "\n")
