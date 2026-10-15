# =============================================================================
# 01. Build per-sample Seurat objects (RNA + ATAC) from Cell Ranger ARC output
#     and compute QC metrics.
# Output: <SEURAT_DIR>/<sample>_raw.rds, QC violin plots
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))

objs <- list()

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]
  h5 <- samples$cnts_h5[i]
  fr <- samples$frags[i]

  cat(glue("[{smp}] Reading 10X h5: {h5}\n"))
  x <- Read10X_h5(h5)

  if (!("Gene Expression" %in% names(x))) {
    stop(glue("[{smp}] Missing 'Gene Expression' in h5. Found: {paste(names(x), collapse=', ')}"))
  }
  if (!("Peaks" %in% names(x))) {
    stop(glue("[{smp}] Missing 'Peaks' in h5. Found: {paste(names(x), collapse=', ')}"))
  }

  cat(glue("[{smp}] Creating SeuratObject... \n"))
  cat(glue("[{smp}] Building RNA assay with Gene Expression... \n"))
  rna_counts <- x[["Gene Expression"]]
  obj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = smp
  )

  cat(glue("[{smp}] Building ChromatinAssay with fragments: {fr}\n"))
  atac_counts <- x[["Peaks"]]
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fr,
    genome = "hg38",
    annotation = annotation
  )
  obj[["ATAC"]] <- chrom_assay

  # Append meta data
  obj$sample <- smp
  obj$condition <- cond

  objs[[smp]] <- obj
  cat(glue("[{smp}] SeuratObject created. \n
            Cells: {ncol(obj)} | Features: RNA - {nrow(obj[['RNA']])}, ATAC - {nrow(obj[['ATAC']])}"))

  cat(glue("[{smp}] Computing basic QC metrics... "))
  obj <- NucleosomeSignal(obj, assay = "ATAC")
  obj <- TSSEnrichment(obj, assay = "ATAC", fast = FALSE)
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")

  cat(glue("[{smp}] Drawing QC metric plots... "))
  suppressMessages({
    library(ggplot2)
    library(cowplot)
    library(gridExtra)
  })

  qc_vals <- c("nFeature_RNA", "nCount_RNA", "percent.mt",
               "nFeature_ATAC", "nCount_ATAC", "nucleosome_signal",
               "TSS.enrichment")

  qc_pannels <- list()

  for (val in qc_vals) {
    # Create QC plot
    qc_pannels[[val]] <- VlnPlot(
      object = obj,
      features = c(val),
      alpha = 0.1,
      pt.size = 0.001
    ) +
      theme(
        plot.title = element_text(size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none"
      )
  }

  qc_plot <- arrangeGrob(
    grobs = qc_pannels,
    ncol = 3
  )

  qc_plot_file <- file.path(qc_plots_dir, paste0(smp, "_qc.png"))
  ggsave(qc_plot_file,
         plot = qc_plot,
         dpi = 600, width = 180, height = 210, unit = "mm")
  cat(glue("[{smp}] QC metric plots saved in: {qc_plot_file}\n"))


  obj_file <- file.path(seurat_dir, paste0(smp, "_raw.rds"))
  cat(glue("[{smp}] Saving SeuratObject as RDS...\n"))
  saveRDS(obj, obj_file)
  cat(glue("[{smp}] Saved in: {obj_file}\n"))
}

cat(glue("[DONE] Raw SeuratObjects and QC plots are created, ready for QC. \n"))
