# =============================================================================
# 01. msVIPER per condition x axis: top 20% vs bottom 20% metacells ranked by
#     axis pressure; gene signature from row-wise t-tests on the full
#     transcriptome CPM matrix (z-scores); msVIPER with minsize = 25.
#
# Input : <RESULTS_DIR>/04_sphingolipid_axes/metacell_obj_with_ceramide_axes_scores_domtype.rds
#         <RESULTS_DIR>/03_viper/mcf7_metacell_regulon_viper.rds
#         <RESULTS_DIR>/01_metacells/metacell_expression_CPM.tsv
# Output: <RESULTS_DIR>/05_axis_regulators/msviper_topbottom20_by_axis/
#           step1_topbottom20_selection.{csv,rds}
#           tb20_<condition>_<axis>_all_genes_{signature.tsv,msviper.rds,results.tsv}
#           all_msviper_topbottom20.rds, all_msviper_topbottom20_results.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(data.table)
  library(Seurat)
  library(SeuratObject)
  library(viper)
  library(genefilter)
})

source(file.path("R", "common.R"))

# Output root
msviper_tb_dir <- step_dir("05_axis_regulators", "msviper_topbottom20_by_axis")

metacell_obj <- readRDS(file.path(DIR_AXES, "metacell_obj_with_ceramide_axes_scores_domtype.rds"))
regulon_obj <- readRDS(file.path(DIR_VIPER, "mcf7_metacell_regulon_viper.rds"))
expr_mat <- read_cpm_matrix()

# Conditions and axes
conditions <- c("control", "acute", "chronic", "return")
axes <- c("S1P", "GSL", "C1P", "SM")

# Top/bottom fraction
top_frac <- 0.20
bottom_frac <- 0.20

# Filters
min_cells_per_group <- 15
minsize_regulon <- 25

# Score columns
axis_score_cols <- c(
  S1P = "Pressure_s1p",
  GSL = "Pressure_gsl",
  C1P = "Pressure_c1p",
  SM  = "Pressure_sm"
)

# Meta
meta <- metacell_obj@meta.data
meta$metacell_id <- colnames(metacell_obj)
meta$condition_std <- tolower(trimws(as.character(meta$condition)))

# Align meta to expr_mat
meta <- meta[match(colnames(expr_mat), meta$metacell_id), , drop = FALSE]
if (!all(meta$metacell_id == colnames(expr_mat))) {
  stop("Metadata order does not match expression matrix columns.")
}

# Check required score columns
need_cols <- c("condition", unname(axis_score_cols))
missing_cols <- setdiff(need_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing meta.data columns: %s", paste(missing_cols, collapse = ", ")))
}

cat("Input check done.\n")
cat("expr_mat:", nrow(expr_mat), "genes x", ncol(expr_mat), "metacells\n")

# -----------------------------------------------------------------------------
# Step 1: top / bottom 20% metacells per condition x axis
# -----------------------------------------------------------------------------
selection_defs <- list()
selection_rows <- list()

for (cond in conditions) {
  meta_cond <- meta[meta$condition_std == cond, , drop = FALSE]

  for (ax in axes) {
    score_col <- axis_score_cols[[ax]]
    x <- meta_cond[[score_col]]

    ok <- is.finite(x)
    meta_sub <- meta_cond[ok, , drop = FALSE]
    x <- x[ok]

    if (length(x) == 0) {
      next
    }

    n_top <- max(1, floor(length(x) * top_frac))
    n_bottom <- max(1, floor(length(x) * bottom_frac))

    ord_desc <- order(x, decreasing = TRUE)
    ord_asc <- order(x, decreasing = FALSE)

    ids_top <- meta_sub$metacell_id[ord_desc[seq_len(n_top)]]
    ids_bottom <- meta_sub$metacell_id[ord_asc[seq_len(n_bottom)]]

    nm <- paste(cond, ax, sep = "__")
    selection_defs[[nm]] <- list(
      condition = cond,
      axis = ax,
      score_col = score_col,
      ids_top = ids_top,
      ids_bottom = ids_bottom
    )

    selection_rows[[nm]] <- data.frame(
      contrast = nm,
      condition = cond,
      axis = ax,
      score_col = score_col,
      n_top = length(ids_top),
      n_bottom = length(ids_bottom),
      top_score_min = min(meta_sub[match(ids_top, meta_sub$metacell_id), score_col], na.rm = TRUE),
      bottom_score_max = max(meta_sub[match(ids_bottom, meta_sub$metacell_id), score_col], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
}

selection_tbl <- do.call(rbind, selection_rows)
write.csv(selection_tbl, file.path(msviper_tb_dir, "step1_topbottom20_selection.csv"), row.names = FALSE)
saveRDS(selection_defs, file.path(msviper_tb_dir, "step1_topbottom20_selection.rds"))

cat("Saved selection table and RDS.\n")

# -----------------------------------------------------------------------------
# Step 2: gene universe (full transcriptome)
# -----------------------------------------------------------------------------
gene_universes <- list(
  all_genes = rownames(expr_mat)
)

# -----------------------------------------------------------------------------
# Step 3: signature (row-wise t-test -> signed z) and msVIPER
# -----------------------------------------------------------------------------
make_signature_topbottom <- function(expr_mat, ids_top, ids_bottom, genes_use) {
  genes_use <- intersect(unique(genes_use), rownames(expr_mat))
  if (length(genes_use) == 0) stop("No overlapping genes in genes_use.")

  x1 <- expr_mat[genes_use, ids_top, drop = FALSE]
  x0 <- expr_mat[genes_use, ids_bottom, drop = FALSE]

  fac <- factor(c(rep("top", ncol(x1)), rep("bottom", ncol(x0))), levels = c("bottom", "top"))
  x <- cbind(x1, x0)

  tt <- genefilter::rowttests(x, fac)

  pval <- as.numeric(tt$p.value)
  stat <- as.numeric(tt$statistic)

  pval[is.na(pval)] <- 1
  pval[pval < 1e-300] <- 1e-300
  pval[pval > 1] <- 1

  zsig <- qnorm(pval / 2, lower.tail = FALSE) * sign(stat)
  names(zsig) <- genes_use

  res <- data.table(
    gene = genes_use,
    statistic = stat,
    p_value = pval,
    z_signature = zsig
  )

  res[, mean_top := rowMeans(x1, na.rm = TRUE)]
  res[, mean_bottom := rowMeans(x0, na.rm = TRUE)]
  res[, diff := mean_top - mean_bottom]
  res[, padj := p.adjust(p_value, method = "BH")]
  res[, abs_z_signature := abs(z_signature)]
  setorder(res, padj, -abs_z_signature)

  list(signature = zsig, stats = res)
}

extract_msviper_results <- function(mra) {
  # summary.msviper() reports the top 10 regulators (default mrs = 10)
  sm <- summary(mra)
  sm_dt <- as.data.table(sm)

  old_names <- colnames(sm_dt)
  new_names <- old_names
  new_names[old_names %in% c("Regulon", "regulon")] <- "regulator"
  new_names[old_names %in% c("Size", "size")] <- "size"
  new_names[old_names %in% c("NES", "nes")] <- "nes"
  new_names[old_names %in% c("p.value", "P.Value", "pvalue")] <- "p.value"
  new_names[old_names %in% c("FDR", "fdr", "padj")] <- "padj"
  setnames(sm_dt, old = old_names, new = new_names)

  if (!("padj" %in% colnames(sm_dt)) && "p.value" %in% colnames(sm_dt)) {
    sm_dt[, padj := p.adjust(p.value, method = "BH")]
  }

  wanted <- c("regulator", "size", "nes", "p.value", "padj")
  keep_first <- wanted[wanted %in% colnames(sm_dt)]
  keep_rest <- setdiff(colnames(sm_dt), keep_first)
  sm_dt[, c(keep_first, keep_rest), with = FALSE]
}

count_regulons_passing_minsize <- function(regulon_obj, genes_use, minsize) {
  genes_use <- unique(genes_use)
  regs <- names(regulon_obj)
  keep <- logical(length(regs))

  for (i in seq_along(regs)) {
    entry <- regulon_obj[[regs[i]]]
    tfm <- entry$tfmode
    if (is.null(tfm)) {
      keep[i] <- FALSE
      next
    }
    ng <- sum(names(tfm) %in% genes_use)
    keep[i] <- (ng >= minsize)
  }

  sum(keep)
}

selection_defs <- readRDS(file.path(msviper_tb_dir, "step1_topbottom20_selection.rds"))

msviper_out_list <- list()
msviper_res_list <- list()

for (nm in names(selection_defs)) {
  sel <- selection_defs[[nm]]
  cond <- sel$condition
  ax <- sel$axis
  ids_top <- sel$ids_top
  ids_bottom <- sel$ids_bottom

  small_n_flag <- (length(ids_top) < min_cells_per_group || length(ids_bottom) < min_cells_per_group)
  if (small_n_flag) {
    message(glue("[RUN small-n] {nm} (top={length(ids_top)}, bottom={length(ids_bottom)})"))
  } else {
    message(glue("[RUN] {nm} (top={length(ids_top)}, bottom={length(ids_bottom)})"))
  }

  universe_map <- list(
    all_genes = gene_universes[["all_genes"]]
  )

  for (univ in names(universe_map)) {
    genes_use <- intersect(unique(universe_map[[univ]]), rownames(expr_mat))
    out_prefix <- paste0("tb20_", cond, "_", ax, "_", univ)

    # Use universe-specific minsize
    minsize_use <- switch(
      univ,
      axis_only = 1L,
      sphingolipid_all = 3L,
      all_genes = minsize_regulon,
      minsize_regulon
    )

    if (length(genes_use) == 0) {
      message(glue("[SKIP no-genes] {out_prefix}"))
      next
    }

    n_reg_ok <- count_regulons_passing_minsize(
      regulon_obj = regulon_obj,
      genes_use = genes_use,
      minsize = minsize_use
    )

    if (n_reg_ok == 0) {
      message(glue("[SKIP no-regulons] {out_prefix} | n_genes={length(genes_use)} | minsize={minsize_use}"))
      next
    }

    message(glue("[RUN universe] {out_prefix} | n_genes={length(genes_use)} | minsize={minsize_use} | regulons={n_reg_ok}"))

    sig_res <- make_signature_topbottom(
      expr_mat = expr_mat,
      ids_top = ids_top,
      ids_bottom = ids_bottom,
      genes_use = genes_use
    )

    fwrite(
      sig_res$stats,
      file.path(msviper_tb_dir, paste0(out_prefix, "_signature.tsv")),
      sep = "\t",
      quote = FALSE
    )

    mra <- tryCatch(
      {
        msviper(
          sig_res$signature,
          regulon_obj,
          minsize = minsize_use,
          ges.filter = TRUE,
          verbose = FALSE
        )
      },
      error = function(e) {
        message(glue("[ERROR msviper] {out_prefix} | {e$message}"))
        NULL
      }
    )

    if (is.null(mra)) next

    saveRDS(mra, file.path(msviper_tb_dir, paste0(out_prefix, "_msviper.rds")))

    res_dt <- extract_msviper_results(mra)
    res_dt[, condition := cond]
    res_dt[, axis := ax]
    res_dt[, universe := univ]
    res_dt[, n_top := length(ids_top)]
    res_dt[, n_bottom := length(ids_bottom)]
    res_dt[, small_n_flag := small_n_flag]
    res_dt[, n_genes_universe := length(genes_use)]
    res_dt[, minsize_used := minsize_use]
    res_dt[, n_regulons_used := n_reg_ok]

    fwrite(
      res_dt,
      file.path(msviper_tb_dir, paste0(out_prefix, "_results.tsv")),
      sep = "\t",
      quote = FALSE
    )

    msviper_out_list[[paste(nm, univ, sep = "__")]] <- mra
    msviper_res_list[[paste(nm, univ, sep = "__")]] <- res_dt
  }
}

all_msviper_tb <- if (length(msviper_res_list) > 0) rbindlist(msviper_res_list, fill = TRUE) else data.table()

saveRDS(msviper_out_list, file.path(msviper_tb_dir, "all_msviper_topbottom20.rds"))
fwrite(
  all_msviper_tb,
  file.path(msviper_tb_dir, "all_msviper_topbottom20_results.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("Saved all msVIPER top-bottom results.\n")
