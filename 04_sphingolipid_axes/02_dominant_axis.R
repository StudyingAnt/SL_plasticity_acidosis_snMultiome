# =============================================================================
# 02. Dominant ceramide clearance axis per metacell (Figure 6c, 6d).
#
#   (1) Threshold diagnostics: 10% / 25% quantiles of
#         range4  = max - min of the four axis pressures
#         margin4 = top1 - top2
#         |p_s1p - p_gsl|
#   (2) Dominant axis = axis with the highest pressure; "Ambiguous" if
#       range4 < range4_thresh or margin4 < margin4_thresh.
#       UMAP colored by softmax-mixed axis colors, stacked bars by condition.
#
# Input : <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores.rds
# Output: <RESULTS_DIR>/04_sphingolipid_axes/threshold_diagnostics/
#         <RESULTS_DIR>/04_sphingolipid_axes/removal_bias_v2_signed/
#           umap_removal_bias_signed_softmax.{svg,png}
#           stack_bias_type_by_condition_signed.{svg,png}
#           stack_dominant4_by_condition_signed.{svg,png}
#         <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores_domtype.rds
# =============================================================================

source(file.path("R", "common.R"))

in_rds <- file.path(DIR_AXES, "metacell_obj_with_ceramide_axes_scores.rds")

if (!file.exists(in_rds)) {
  stop(paste("File not found:", in_rds))
}

metacell_obj <- readRDS(in_rds)

# -----------------------------------------------------------------------------
# (1) Threshold diagnostics
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(SeuratObject)
})

# Output directory
sub_dir <- step_dir("04_sphingolipid_axes", "threshold_diagnostics")

# Pull meta.data
md <- metacell_obj[[]]
need_cols <- c("Pressure_s1p", "Pressure_gsl", "Pressure_c1p", "Pressure_sm")
missing_cols <- setdiff(need_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing meta.data columns: %s", paste(missing_cols, collapse = ", ")))
}

df <- data.frame(
  p_s1p = md$Pressure_s1p,
  p_gsl = md$Pressure_gsl,
  p_c1p = md$Pressure_c1p,
  p_sm  = md$Pressure_sm,
  stringsAsFactors = FALSE
)

# Compute diagnostics
pmax4 <- pmax(df$p_s1p, pmax(df$p_gsl, pmax(df$p_c1p, df$p_sm)))
pmin4 <- pmin(df$p_s1p, pmin(df$p_gsl, pmin(df$p_c1p, df$p_sm)))
df$range4 <- pmax4 - pmin4

P4 <- as.matrix(df[, c("p_s1p","p_gsl","p_c1p","p_sm"), drop = FALSE])
top1_idx <- max.col(P4, ties.method = "first")
top1 <- P4[cbind(seq_len(nrow(P4)), top1_idx)]
top2 <- apply(P4, 1, function(v) sort(v, decreasing = TRUE)[2])
df$margin4 <- top1 - top2

df$diff_s1p_gsl <- df$p_s1p - df$p_gsl
df$abs_diff_s1p_gsl <- abs(df$diff_s1p_gsl)

# Drop non-finite
df <- df[is.finite(df$range4) & is.finite(df$margin4) & is.finite(df$abs_diff_s1p_gsl), , drop = FALSE]

# Quantile thresholds to report (10% and 25%)
q_levels <- c(0.10, 0.25)

thr_range4 <- quantile(df$range4, probs = q_levels, na.rm = TRUE, names = TRUE)
thr_margin4 <- quantile(df$margin4, probs = q_levels, na.rm = TRUE, names = TRUE)
thr_absdiff <- quantile(df$abs_diff_s1p_gsl, probs = q_levels, na.rm = TRUE, names = TRUE)

# Print thresholds
cat("\n=== Quantile thresholds (10% and 25%) ===\n")
cat("range4 (max-min):\n"); print(thr_range4)
cat("\nmargin4 (top1-top2):\n"); print(thr_margin4)
cat("\nabs(diff_s1p_gsl):\n"); print(thr_absdiff)

# Save thresholds to a text file
thr_file <- file.path(sub_dir, "thresholds_q10_q25.txt")
writeLines(
  c(
    "Quantile thresholds (10% and 25%)",
    "",
    "range4 = max(p) - min(p)",
    paste(names(thr_range4), sprintf("%.6f", as.numeric(thr_range4)), sep = "\t"),
    "",
    "margin4 = top1(p) - top2(p)",
    paste(names(thr_margin4), sprintf("%.6f", as.numeric(thr_margin4)), sep = "\t"),
    "",
    "abs_diff_s1p_gsl = |p_s1p - p_gsl|",
    paste(names(thr_absdiff), sprintf("%.6f", as.numeric(thr_absdiff)), sep = "\t")
  ),
  con = thr_file
)
cat("\nSaved threshold table to:", thr_file, "\n")

# Helper: plot distribution + vertical lines for q10/q25
plot_dist_with_thresholds <- function(x, thr_named, title, xlab, bins = 60) {
  d <- data.frame(x = x)
  v10 <- as.numeric(thr_named[1])
  v25 <- as.numeric(thr_named[2])

  ggplot(d, aes(x = x)) +
    geom_histogram(bins = bins, color = "black", linewidth = 0.2) +
    geom_vline(xintercept = v10, linewidth = 0.45, linetype = "dashed") +
    geom_vline(xintercept = v25, linewidth = 0.45, linetype = "dotted") +
    labs(
      title = title,
      x = xlab,
      y = "Count"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
}

# Make plots
p_range4 <- plot_dist_with_thresholds(
  x = df$range4,
  thr_named = thr_range4,
  title = "range4 = max(p) - min(p) (all cells)",
  xlab = "range4",
  bins = 70
)

p_margin4 <- plot_dist_with_thresholds(
  x = df$margin4,
  thr_named = thr_margin4,
  title = "margin4 = top1(p) - top2(p) (all cells)",
  xlab = "margin4",
  bins = 70
)

p_absdiff <- plot_dist_with_thresholds(
  x = df$abs_diff_s1p_gsl,
  thr_named = thr_absdiff,
  title = "|diff_s1p_gsl| = |p_s1p - p_gsl| (all cells)",
  xlab = "|p_s1p - p_gsl|",
  bins = 70
)

# Save plots
ggsave(file.path(sub_dir, "dist_range4_with_q10_q25.svg"),
       p_range4, width = 3.8, height = 2.6, units = "in", bg = "white")
ggsave(file.path(sub_dir, "dist_range4_with_q10_q25.png"),
       p_range4, width = 3.8, height = 2.6, units = "in", dpi = 600, bg = "white")

ggsave(file.path(sub_dir, "dist_margin4_with_q10_q25.svg"),
       p_margin4, width = 3.8, height = 2.6, units = "in", bg = "white")
ggsave(file.path(sub_dir, "dist_margin4_with_q10_q25.png"),
       p_margin4, width = 3.8, height = 2.6, units = "in", dpi = 600, bg = "white")

ggsave(file.path(sub_dir, "dist_absdiff_s1p_gsl_with_q10_q25.svg"),
       p_absdiff, width = 3.8, height = 2.6, units = "in", bg = "white")
ggsave(file.path(sub_dir, "dist_absdiff_s1p_gsl_with_q10_q25.png"),
       p_absdiff, width = 3.8, height = 2.6, units = "in", dpi = 600, bg = "white")

cat("Saved distribution plots to:", sub_dir, "\n")

# -----------------------------------------------------------------------------
# (2) Dominant axis, UMAP and stacked bars
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(grid)   # <-- needed for unit()
})

# =========================
# Output / plotting settings
# =========================
sub_dir <- step_dir("04_sphingolipid_axes", "removal_bias_v2_signed")

reduction_name <- "umap_sct_mc"
pt_size <- 0.25
pt_alpha <- 0.95

# ---- minimal safety: reduction fallback (only if missing)
if (!reduction_name %in% names(metacell_obj@reductions)) {
  umap_candidates <- grep("^umap", names(metacell_obj@reductions), value = TRUE)
  if (length(umap_candidates) > 0) {
    reduction_name <- umap_candidates[1]
  } else {
    stop("No UMAP reduction found in metacell_obj@reductions.")
  }
}

# =========================
# Thresholds
# =========================
# Values from threshold_diagnostics/thresholds_q10_q25.txt (10% quantiles)
diff_switch <- 0.1862378
range4_thresh <- 1.793752
margin4_thresh <- 0.1855446
beta <- 1.0

# =========================
# Colors
# =========================
col_s1p <- "#CC79A7"
col_gsl <- "#009E73"
col_c1p <- "#E69F00"
col_sm  <- "#0072B2"
col_ambig <- "#999999"

hex_to_rgb01 <- function(hex) {
  x <- grDevices::col2rgb(hex)
  as.numeric(x) / 255
}

rgb_s1p <- hex_to_rgb01(col_s1p)
rgb_gsl <- hex_to_rgb01(col_gsl)
rgb_c1p <- hex_to_rgb01(col_c1p)
rgb_sm  <- hex_to_rgb01(col_sm)

softmax3 <- function(a, b, c, beta = 1.0) {
  m <- pmax(a, pmax(b, c))
  ea <- exp(beta * (a - m))
  eb <- exp(beta * (b - m))
  ec <- exp(beta * (c - m))
  s <- ea + eb + ec
  cbind(ea / s, eb / s, ec / s)
}

mix3_softmax <- function(x1, x2, x3, col1, col2, col3, beta = 1.0) {
  W <- softmax3(x1, x2, x3, beta = beta)
  r <- W[, 1] * col1[1] + W[, 2] * col2[1] + W[, 3] * col3[1]
  g <- W[, 1] * col1[2] + W[, 2] * col2[2] + W[, 3] * col3[2]
  b <- W[, 1] * col1[3] + W[, 2] * col2[3] + W[, 3] * col3[3]
  grDevices::rgb(r, g, b)
}

# =========================
# Pull data
# =========================
metacell_obj$condition <- factor(
  metacell_obj$condition,
  levels = CONDITION_LEVELS
)

md <- metacell_obj[[]]
need_cols <- c("condition", "Pressure_s1p", "Pressure_gsl", "Pressure_c1p", "Pressure_sm")
missing_cols <- setdiff(need_cols, colnames(md))
if (length(missing_cols) > 0) stop(sprintf("Missing meta.data columns: %s", paste(missing_cols, collapse = ", ")))

um <- Embeddings(metacell_obj, reduction = reduction_name)
umdf <- data.frame(cell = rownames(um), UMAP_1 = um[, 1], UMAP_2 = um[, 2], stringsAsFactors = FALSE)

df <- data.frame(
  cell = rownames(md),
  condition = md$condition,
  p_s1p = md$Pressure_s1p,
  p_gsl = md$Pressure_gsl,
  p_c1p = md$Pressure_c1p,
  p_sm  = md$Pressure_sm,
  stringsAsFactors = FALSE
)

df <- merge(df, umdf, by = "cell", all.x = TRUE)
df <- df[!is.na(df$UMAP_1) & !is.na(df$UMAP_2), , drop = FALSE]

# =========================
# Fully ambiguous definition from signed pressures
# =========================
pmax4 <- pmax(df$p_s1p, pmax(df$p_gsl, pmax(df$p_c1p, df$p_sm)))
pmin4 <- pmin(df$p_s1p, pmin(df$p_gsl, pmin(df$p_c1p, df$p_sm)))
df$range4 <- pmax4 - pmin4

# =========================
# Build two view colors using SOFTMAX on signed pressures
# =========================
df$col_view_s1p <- mix3_softmax(df$p_s1p, df$p_c1p, df$p_sm, rgb_s1p, rgb_c1p, rgb_sm, beta = beta)
df$col_view_gsl <- mix3_softmax(df$p_gsl, df$p_c1p, df$p_sm, rgb_gsl, rgb_c1p, rgb_sm, beta = beta)

# =========================
# bias_type: decide which view to use, plus ambiguous handling
# =========================
df$diff_s1p_gsl <- df$p_s1p - df$p_gsl

df$bias_type <- NA_character_
df$bias_type[df$range4 < range4_thresh] <- "FullyAmbig"

idx_ok <- which(df$range4 >= range4_thresh)
if (length(idx_ok) > 0) {
  d <- df$diff_s1p_gsl[idx_ok]
  df$bias_type[idx_ok[d >=  diff_switch]] <- "ViewS1P"
  df$bias_type[idx_ok[d <= -diff_switch]] <- "ViewGSL"

  idx_mid <- idx_ok[abs(d) < diff_switch]
  if (length(idx_mid) > 0) {
    df$bias_type[idx_mid[df$p_s1p[idx_mid] >= df$p_gsl[idx_mid]]] <- "Ambig_S1P"
    df$bias_type[idx_mid[df$p_s1p[idx_mid] <  df$p_gsl[idx_mid]]] <- "Ambig_GSL"
  }
}

df$bias_type <- factor(
  df$bias_type,
  levels = c("ViewS1P", "ViewGSL", "Ambig_S1P", "Ambig_GSL", "FullyAmbig")
)

df$color_final <- NA_character_
df$color_final[df$bias_type == "ViewS1P"] <- df$col_view_s1p[df$bias_type == "ViewS1P"]
df$color_final[df$bias_type == "ViewGSL"] <- df$col_view_gsl[df$bias_type == "ViewGSL"]
df$color_final[df$bias_type == "Ambig_S1P"] <- df$col_view_s1p[df$bias_type == "Ambig_S1P"]
df$color_final[df$bias_type == "Ambig_GSL"] <- df$col_view_gsl[df$bias_type == "Ambig_GSL"]
df$color_final[df$bias_type == "FullyAmbig"] <- col_ambig

# =========================
# UMAP plot
# =========================
p_umap <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = df$color_final, size = pt_size, alpha = pt_alpha) +
  coord_equal() +
  labs(
    title = sprintf("Removal bias (signed pressure softmax mix; beta=%.2f)", beta),
    x = "UMAP_1",
    y = "UMAP_2"
  ) +
  theme_classic(base_size = 7) +
  theme(
    plot.title = element_text(size = 7, hjust = 0.5),
    axis.title = element_text(size = 7),
    axis.text  = element_text(size = 6),
    aspect.ratio = 1
  )

save_plot_pair(p_umap, "umap_removal_bias_signed_softmax", sub_dir, w = 2.0, h = 2.0, dpi = 600)

# =========================
# Stacked bar 1: bias_type proportion by condition (100%)
# =========================
type_colors <- c(
  "ViewS1P"    = col_s1p,
  "ViewGSL"    = col_gsl,
  "Ambig_S1P"  = "#9B59B6",
  "Ambig_GSL"  = "#4CAF50",
  "FullyAmbig" = col_ambig
)

p_stack1 <- ggplot(df, aes(x = condition, fill = bias_type)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = type_colors, drop = FALSE) +
  scale_y_continuous(
    name = "Proportion",
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = sprintf("Bias-type proportion | diff_switch=%.2f | range4<%.2f => FullyAmbig", diff_switch, range4_thresh),
    x = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  )

save_plot_pair(p_stack1, "stack_bias_type_by_condition_signed", sub_dir, w = 2.2, h = 1.3, dpi = 600)

# =========================
# Stacked bar 2: dominant among 4 signed pressures + ambiguous
# =========================
P4 <- as.matrix(df[, c("p_s1p","p_gsl","p_c1p","p_sm"), drop = FALSE])

top1_idx <- max.col(P4, ties.method = "first")
top1 <- P4[cbind(seq_len(nrow(P4)), top1_idx)]
top2 <- apply(P4, 1, function(v) sort(v, decreasing = TRUE)[2])
margin4 <- top1 - top2

dominant <- c("S1P","GSL","C1P","SM")[top1_idx]
dominant[df$range4 < range4_thresh] <- "Ambiguous"
dominant[margin4 < margin4_thresh] <- "Ambiguous"

df$dominant4 <- factor(dominant, levels = c("S1P","GSL","C1P","SM","Ambiguous"))

dom_colors <- c(
  "S1P"       = col_s1p,
  "GSL"       = col_gsl,
  "C1P"       = col_c1p,
  "SM"        = col_sm,
  "Ambiguous" = col_ambig
)

p_stack2 <- ggplot(df, aes(x = condition, fill = dominant4)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = dom_colors, drop = FALSE) +
  scale_y_continuous(
    name = "Proportion",
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) x * 100,
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = sprintf("Dominant axis among 4 (signed pressure) | margin4<%.2f or range4<%.2f => Ambiguous", margin4_thresh, range4_thresh),
    x = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  )

save_plot_pair(p_stack2, "stack_dominant4_by_condition_signed", sub_dir, w = 3.0, h = 1.8, dpi = 600)

cat("Saved signed-pressure UMAP + stacked bars to:", sub_dir, "\n")

# -----------------------------------------------------------------------------
# Store dominant axis in the object
# -----------------------------------------------------------------------------
# df was re-ordered by merge() above, so labels are matched to metacells by name.
# (The original notebook assigned the vector by position,
#  metacell_obj$dominant4 <- dominant. No figure in this repository uses the
#  stored dominant4 column; Figure 6c/6d are drawn from df above.)
metacell_obj$dominant4 <- setNames(dominant, df$cell)[colnames(metacell_obj)]

out_rds <- file.path(DIR_AXES, "metacell_obj_with_ceramide_axes_scores_domtype.rds")
saveRDS(metacell_obj, file = out_rds)
cat("Saved updated object to:", out_rds, "\n")
