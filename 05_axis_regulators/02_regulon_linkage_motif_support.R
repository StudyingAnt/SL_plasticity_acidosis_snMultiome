# =============================================================================
# 02. Integrate msVIPER results with
#       (1) axis gene linkage: sum |tfmode * likelihood| from each regulator
#           to the axis genes in the regulon, and
#       (2) motif accessibility support: Wilcoxon test of metacell chromVAR
#           deviations (top vs bottom 20%) for JASPAR motifs of each TF.
#
# Input : <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores_domtype.rds
#         <RESULTS_DIR>/03_viper/mcf7_metacell_regulon_viper.rds
#         <RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/ (step 01)
# Output: <RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/
#           regulon_axis_stats.tsv
#           all_msviper_topbottom20_with_regulon.tsv
#           motif_validation_{summary,detail}.tsv
#           final_candidates_topbottom20.tsv (+ per condition x axis tables)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(SeuratObject)
})

source(file.path("R", "common.R"))

msviper_tb_dir <- file.path(DIR_REGS, "msviper_topbottom20_by_axis")

metacell_obj <- readRDS(file.path(DIR_AXES, "metacell_obj_with_ceramide_axes_scores_domtype.rds"))
regulon_obj <- readRDS(file.path(DIR_VIPER, "mcf7_metacell_regulon_viper.rds"))

axes <- c("S1P", "GSL", "C1P", "SM")

# Axis gene sets (R/common.R)
axis_gene_sets <- list(
  S1P = s1p_fin_axis,
  GSL = gsl_axis,
  C1P = c1p_axis,
  SM  = sm_axis
)

# Number of top msVIPER regulators per contrast passed to motif validation
top_n_tf_for_motif <- 50

# -----------------------------------------------------------------------------
# Regulon linkage to axis genes
# -----------------------------------------------------------------------------
regulon_axis_stats <- function(regulon_obj, axis_genes) {
  regs <- names(regulon_obj)
  axis_genes <- unique(axis_genes)

  out <- data.table(
    regulator = regs,
    n_targets_in_axis = 0L,
    sum_abs_w = 0,
    sum_w = 0,
    max_abs_w = 0
  )

  for (i in seq_along(regs)) {
    r <- regs[i]
    entry <- regulon_obj[[r]]
    tfm <- entry$tfmode
    lik <- entry$likelihood
    if (is.null(tfm) || is.null(lik)) next

    idx <- match(axis_genes, names(tfm))
    ok <- which(!is.na(idx))
    if (length(ok) == 0) next

    w <- unname(tfm[idx[ok]]) * unname(lik[idx[ok]])
    out$n_targets_in_axis[i] <- length(w)
    out$sum_abs_w[i] <- sum(abs(w))
    out$sum_w[i] <- sum(w)
    out$max_abs_w[i] <- max(abs(w))
  }

  out
}

regulon_axis_all <- rbindlist(lapply(axes, function(ax) {
  tab <- regulon_axis_stats(regulon_obj, axis_gene_sets[[ax]])
  tab[, axis := ax]
  tab
}), fill = TRUE)

fwrite(regulon_axis_all, file.path(msviper_tb_dir, "regulon_axis_stats.tsv"), sep = "\t", quote = FALSE)

tb <- fread(file.path(msviper_tb_dir, "all_msviper_topbottom20_results.tsv"))
tb[, axis := toupper(axis)]

tb_reg <- merge(
  tb,
  regulon_axis_all,
  by = c("axis", "regulator"),
  all.x = TRUE
)

tb_reg[is.na(n_targets_in_axis), n_targets_in_axis := 0L]
tb_reg[is.na(sum_abs_w), sum_abs_w := 0]
tb_reg[is.na(sum_w), sum_w := 0]
tb_reg[is.na(max_abs_w), max_abs_w := 0]

tb_reg[, abs_nes := abs(nes)]
tb_reg[, net_strength := log1p(sum_abs_w)]
tb_reg[, net_targets := log1p(n_targets_in_axis)]
tb_reg[, net_dir := fifelse(sum_w > 0, 1, fifelse(sum_w < 0, -1, 0))]
tb_reg[, aligned_nes := nes * net_dir]

tb_reg[, score_net_nes :=
  1.0 * abs_nes +
  0.6 * net_strength +
  0.3 * net_targets +
  0.5 * aligned_nes
]

setorder(tb_reg, universe, condition, axis, -score_net_nes, padj)

fwrite(tb_reg, file.path(msviper_tb_dir, "all_msviper_topbottom20_with_regulon.tsv"),
       sep = "\t", quote = FALSE)

cat("Saved msVIPER + regulon table.\n")

# -----------------------------------------------------------------------------
# Motif accessibility support (metacell chromVAR)
# -----------------------------------------------------------------------------
motif_name_map <- metacell_obj@misc$motif.names
chromvar_mat <- as.matrix(GetAssayData(metacell_obj, assay = "chromvar_mc", layer = "data"))

split_tf_tokens <- function(motif_name) {
  x <- toupper(motif_name)
  x <- gsub("\\s+", "", x)
  parts <- unlist(strsplit(x, "::", fixed = TRUE))
  parts <- unlist(strsplit(parts, "/", fixed = TRUE))
  parts <- unlist(strsplit(parts, ":", fixed = TRUE))
  parts <- parts[parts != ""]
  unique(parts)
}

motifs_for_tf <- function(tf, motif_name_map) {
  tfu <- toupper(tf)
  nm <- as.character(motif_name_map)
  hit <- logical(length(nm))
  for (i in seq_along(nm)) {
    toks <- split_tf_tokens(nm[i])
    if (tfu %in% toks) hit[i] <- TRUE
  }
  names(motif_name_map)[hit]
}

chromvar_de_one <- function(mat, cells.1, cells.2, features) {
  c1 <- intersect(cells.1, colnames(mat))
  c2 <- intersect(cells.2, colnames(mat))
  feats <- intersect(features, rownames(mat))

  if (length(c1) == 0 || length(c2) == 0 || length(feats) == 0) return(NULL)

  x1 <- mat[feats, c1, drop = FALSE]
  x0 <- mat[feats, c2, drop = FALSE]

  pval <- vapply(seq_len(nrow(x1)), function(i) {
    a <- x1[i, ]
    b <- x0[i, ]
    a <- a[is.finite(a)]
    b <- b[is.finite(b)]
    if (length(unique(c(a, b))) <= 1) return(1)
    suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  }, numeric(1))

  mean1 <- rowMeans(x1, na.rm = TRUE)
  mean0 <- rowMeans(x0, na.rm = TRUE)

  out <- data.table(
    motif_id = feats,
    avg_diff = mean1 - mean0,
    p_val = pval,
    p_val_adj = p.adjust(pval, method = "BH")
  )
  out
}

selection_defs <- readRDS(file.path(msviper_tb_dir, "step1_topbottom20_selection.rds"))
tb_reg <- fread(file.path(msviper_tb_dir, "all_msviper_topbottom20_with_regulon.tsv"))

motif_summary_list <- list()
motif_detail_list <- list()

for (nm in names(selection_defs)) {
  sel <- selection_defs[[nm]]
  cond <- sel$condition
  ax <- sel$axis
  ids_top <- sel$ids_top
  ids_bottom <- sel$ids_bottom

  for (univ in c("all_genes")) {
    cand <- tb_reg[
      condition == cond & axis == ax & universe == univ
    ][order(-score_net_nes)][1:min(.N, top_n_tf_for_motif)]

    if (nrow(cand) == 0) next

    for (tf in cand$regulator) {
      mids <- motifs_for_tf(tf, motif_name_map)
      if (length(mids) == 0) next

      de <- chromvar_de_one(
        mat = chromvar_mat,
        cells.1 = ids_top,
        cells.2 = ids_bottom,
        features = mids
      )
      if (is.null(de) || nrow(de) == 0) next

      de[, regulator := tf]
      de[, condition := cond]
      de[, axis := ax]
      de[, universe := univ]
      motif_detail_list[[paste(cond, ax, univ, tf, sep = "__")]] <- de

      motif_summary_list[[paste(cond, ax, univ, tf, sep = "__")]] <- data.table(
        regulator = tf,
        condition = cond,
        axis = ax,
        universe = univ,
        motif_n = nrow(de),
        motif_best_padj = min(de$p_val_adj, na.rm = TRUE),
        motif_max_abs_diff = max(abs(de$avg_diff), na.rm = TRUE)
      )
    }
  }
}

motif_summary_dt <- if (length(motif_summary_list) > 0) rbindlist(motif_summary_list, fill = TRUE) else data.table()
motif_detail_dt <- if (length(motif_detail_list) > 0) rbindlist(motif_detail_list, fill = TRUE) else data.table()

fwrite(motif_summary_dt, file.path(msviper_tb_dir, "motif_validation_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(motif_detail_dt, file.path(msviper_tb_dir, "motif_validation_detail.tsv"), sep = "\t", quote = FALSE)

cat("Saved motif validation summary/detail.\n")

# -----------------------------------------------------------------------------
# Final candidate table
# -----------------------------------------------------------------------------
tb_reg <- fread(file.path(msviper_tb_dir, "all_msviper_topbottom20_with_regulon.tsv"))
motif_summary_dt <- fread(file.path(msviper_tb_dir, "motif_validation_summary.tsv"))

final_dt <- merge(
  tb_reg,
  motif_summary_dt,
  by = c("condition", "axis", "universe", "regulator"),
  all.x = TRUE
)

final_dt[is.na(motif_n), motif_n := 0L]
final_dt[is.na(motif_best_padj), motif_best_padj := 1]
final_dt[is.na(motif_max_abs_diff), motif_max_abs_diff := 0]

final_dt[, motif_score := -log10(motif_best_padj + 1e-300) + 0.25 * motif_max_abs_diff]

final_dt[, final_score :=
  score_net_nes + motif_score
]

setorder(final_dt, universe, condition, axis, -final_score, padj)

fwrite(final_dt, file.path(msviper_tb_dir, "final_candidates_topbottom20.tsv"),
       sep = "\t", quote = FALSE)

for (univ in unique(final_dt$universe)) {
  for (cond in unique(final_dt$condition)) {
    for (ax in unique(final_dt$axis)) {
      sub <- final_dt[universe == univ & condition == cond & axis == ax]
      if (nrow(sub) == 0) next
      fwrite(
        sub,
        file.path(msviper_tb_dir, paste0("final_", univ, "_", cond, "_", ax, ".tsv")),
        sep = "\t",
        quote = FALSE
      )
    }
  }
}

cat("Saved final candidate tables.\n")
