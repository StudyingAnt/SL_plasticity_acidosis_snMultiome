# =============================================================================
# 01. Extended Data figure: where the ATAC data enter the analysis.
#
#   a-c  UMAP of all nuclei colored by condition, computed from
#        a: RNA (SCT PCA), b: ATAC (peaks LSI), c: WNN (RNA + ATAC)
#   d    Mean chromVAR deviation per condition for the JASPAR motifs of the
#        candidate regulators shown in Figure 6f-i (same selection as the dot
#        plots: msVIPER padj < 0.05, top 20 per axis by mean -log10 padj)
#
# Input : <SEURAT_DIR>/all_qc_noMT_macs3_trf_chrvar.rds (00-2 step 07)
#         <RESULTS_DIR>/05_axis_regulators/.../plot_dt_allgenes_regulon_motifalpha_padj005.tsv
# Output: <RESULTS_DIR>/06_extended_data_atac/
#           umap_{rna,atac,wnn}_by_condition.{svg,png}, umap_legend.{svg,png}
#           chromvar_candidate_regulators_heatmap.{svg,png}
#           chromvar_candidate_regulators_mean_by_condition.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Signac)
  library(data.table)
  library(ggplot2)
})

source(file.path("R", "common.R"))

out_dir <- step_dir("06_extended_data_atac")

obj <- readRDS(file.path(SEURAT_DIR, "all_qc_noMT_macs3_trf_chrvar.rds"))
obj$condition <- factor(obj$condition, levels = CONDITION_LEVELS)
cat("Nuclei:", ncol(obj), "| reductions:", paste(Reductions(obj), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# a-c: UMAPs by condition
# -----------------------------------------------------------------------------
umaps <- c(rna = "rna.sct.umap", atac = "peaks.umap", wnn = "mm.sct.peaks.umap")
titles <- c(rna = "RNA (SCT)", atac = "ATAC (peaks LSI)", wnn = "WNN (RNA + ATAC)")

for (nm in names(umaps)) {
  # raster: pt.size is in pixels of the raster.dpi grid; shuffle so no condition is drawn on top
  p <- DimPlot(obj, reduction = umaps[[nm]], group.by = "condition", shuffle = TRUE, seed = 1234,
               pt.size = 1.6, raster = TRUE, raster.dpi = c(800, 800)) +
    scale_color_manual(values = CONDITION_COLORS, breaks = CONDITION_LEVELS, drop = FALSE) +
    labs(title = titles[[nm]], x = "UMAP_1", y = "UMAP_2", color = "Condition") +
    theme_classic(base_size = 7) +
    theme(plot.title = element_text(size = 7, hjust = 0.5),
          axis.title = element_text(size = 7), axis.text = element_text(size = 6),
          legend.position = "none", aspect.ratio = 1)
  save_plot_pair(p, paste0("umap_", nm, "_by_condition"), out_dir, w = 1.8, h = 1.9)
  if (nm == "rna") {
    save_plot_pair(p + theme(legend.position = "right", legend.title = element_text(size = 7),
                             legend.text = element_text(size = 6)),
                   "umap_legend", out_dir, w = 2.6, h = 1.9)
  }
}

# -----------------------------------------------------------------------------
# d: chromVAR deviations of candidate regulators (Figure 6f-i selection)
# -----------------------------------------------------------------------------
plot_file <- file.path(DIR_REGS, "msviper_topbottom20_by_axis", "viz_allgenes_regulon_motifalpha_padj005",
                       "plot_dt_allgenes_regulon_motifalpha_padj005.tsv")
pd <- fread(plot_file)
pd[, axis := toupper(axis)]
pd <- pd[padj < 0.05]
topM <- 20
sel <- rbindlist(lapply(unique(pd$axis), function(ax) {
  d <- pd[axis == ax, .(mean_msviper_logp = mean(-log10(padj + 1e-300), na.rm = TRUE)), by = regulator]
  d <- d[order(-mean_msviper_logp)][seq_len(min(.N, topM))]
  d[, axis := ax]
}))

# same TF -> JASPAR motif matching as 05_axis_regulators/02_regulon_linkage_motif_support.R
split_tf_tokens <- function(motif_name) {
  x <- toupper(motif_name)
  x <- gsub("\\s+", "", x)
  parts <- unlist(strsplit(x, "::", fixed = TRUE))
  parts <- unlist(strsplit(parts, "/", fixed = TRUE))
  parts <- unlist(strsplit(parts, ":", fixed = TRUE))
  parts <- parts[parts != ""]
  unique(parts)
}
motif_names <- obj@assays[["peaks"]]@motifs@motif.names
motifs_for_tf <- function(tf) {
  hit <- vapply(as.character(motif_names), function(n) toupper(tf) %in% split_tf_tokens(n), logical(1))
  names(motif_names)[hit]
}

cv <- GetAssayData(obj, assay = "chromvar", layer = "data")
cond <- obj$condition

rows <- list()
for (i in seq_len(nrow(sel))) {
  tf <- sel$regulator[i]
  for (m in intersect(motifs_for_tf(tf), rownames(cv))) {
    z <- tapply(as.numeric(cv[m, ]), cond, mean, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.table(axis = sel$axis[i], regulator = tf, motif_id = m,
                                           motif_name = motif_names[[m]],
                                           condition = names(z), mean_deviation = as.numeric(z))
  }
}
hm <- rbindlist(rows)
fwrite(hm, file.path(out_dir, "chromvar_candidate_regulators_mean_by_condition.tsv"), sep = "\t")

n_tf <- sel[, .N, by = axis]
n_motif <- unique(hm[, .(axis, regulator)])[, .N, by = axis]
cat("Candidate regulators per axis:", paste(n_tf$axis, n_tf$N, collapse = ", "), "\n")
cat("  ... with a JASPAR motif:", paste(n_motif$axis, n_motif$N, collapse = ", "), "\n")

hm[, label := paste0(regulator, " (", motif_id, ")")]
hm[, condition := factor(condition, levels = CONDITION_LEVELS)]
hm[, axis := factor(axis, levels = c("S1P", "GSL", "C1P", "SM"))]
# y key is unique per axis, so a regulator shared by several axes is sorted within each panel
hm[, row_key := paste(axis, label, sep = "|")]
ord <- unique(hm[order(axis, regulator, motif_id), row_key])
hm[, row_key := factor(row_key, levels = rev(ord))]
lim <- max(abs(hm$mean_deviation), na.rm = TRUE)

p_hm <- ggplot(hm, aes(x = condition, y = row_key, fill = mean_deviation)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "Mean chromVAR\ndeviation (z)") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 5.5),
        strip.text.y = element_text(size = 7, angle = 0),
        legend.title = element_text(size = 6), legend.text = element_text(size = 6))

h <- max(2.5, 0.09 * length(unique(hm$row_key)) + 0.8)
save_plot_pair(p_hm, "chromvar_candidate_regulators_heatmap", out_dir, w = 3.0, h = h)
cat("Saved figures to", out_dir, "\n")
