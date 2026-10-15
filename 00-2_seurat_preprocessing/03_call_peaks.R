# =============================================================================
# 03. Call peaks per sample with MACS3, remove non-standard chromosomes and
#     hg38 blacklist regions, build a union peak set, and add a unified
#     'peaks' ChromatinAssay to each sample. Also exports the peak matrices.
# Output: <SEURAT_DIR>/<sample>_qc_noMT_macs3.rds
# =============================================================================

source(file.path("00-2_seurat_preprocessing", "00_setup.R"))
suppressMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# MACS3 binary (MACS3_PATH) is set in config/config.R
# macs_out_root is set in 00_setup.R
# blacklist_hg38_unified is provided by Signac

cat(glue("MACS3_PATH: {MACS3_PATH}\n"))
cat(glue("macs_out_root: {macs_out_root}\n\n"))

# Create MACS3 output root
dir.create(macs_out_root, recursive = TRUE, showWarnings = FALSE)

objs <- list()
peaks_ls <- list()

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]

  obj_file <- file.path(seurat_dir, paste0(smp, "_qc_noMT.rds"))
  cat(glue("[{smp}] Loading SeuratObject from: {obj_file}\n"))
  obj <- readRDS(obj_file)

  cat(glue(
    "[{smp}] SeuratObject loaded.\n",
    "  Condition: {cond}\n",
    "  Assays: {paste(Assays(obj), collapse = ', ')}\n",
    "  Cells: {ncol(obj)}\n",
    "  Fragments: {length(Fragments(obj[['ATAC']]))}\n"
  ))

  objs[[smp]] <- obj

  outdir_sid <- file.path(macs_out_root, smp)
  dir.create(outdir_sid, recursive = TRUE, showWarnings = FALSE)

  cat(glue("[{smp}] Calling peaks with MACS3...\n"))
  cat(glue("[{smp}] MACS3 output directory: {outdir_sid}\n"))

  peaks <- CallPeaks(
    object = obj,
    assay = "ATAC",
    macs2.path = MACS3_PATH,
    outdir = outdir_sid,
    cleanup = FALSE,
    verbose = TRUE
  )

  cat(glue(
    "[{smp}] Peak calling completed.\n",
    "  Raw peaks: {length(peaks)}\n"
  ))

  cat(glue("[{smp}] Filtering peaks to standard chromosomes...\n"))
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  cat(glue("[{smp}] Peaks after standard chromosome filter: {length(peaks)}\n"))

  cat(glue("[{smp}] Removing blacklist regions...\n"))
  peaks <- subsetByOverlaps(
    x = peaks,
    ranges = blacklist_hg38_unified,
    invert = TRUE
  )
  cat(glue("[{smp}] Peaks after blacklist removal: {length(peaks)}\n"))

  peaks_ls[[smp]] <- peaks

  peaks_rds <- file.path(
    macs_out_root,
    paste0(smp, "_macs3_peaks.rds")
  )
  saveRDS(peaks, peaks_rds)
  cat(glue("[{smp}] Saved MACS3 peaks in: {peaks_rds}\n\n"))
}

cat(glue("[UNION] Creating union peaks across all samples...\n"))

all_peaks <- unlist(GenomicRanges::GRangesList(peaks_ls), use.names = FALSE)
peaks_union <- GenomicRanges::reduce(all_peaks, ignore.strand = TRUE)

cat(glue(
  "[UNION] Union peaks created.\n",
  "  Number of union peaks: {length(peaks_union)}\n"
))

peaks_union_rds <- file.path(macs_out_root, "peaks_union.rds")
saveRDS(peaks_union, peaks_union_rds)
cat(glue("[UNION] Saved union peaks in: {peaks_union_rds}\n\n"))

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]

  cat(glue("[{smp}] Building unified peaks assay...\n"))
  obj <- objs[[smp]]
  DefaultAssay(obj) <- "ATAC"

  mat <- FeatureMatrix(
    fragments = Fragments(obj),
    features = peaks_union,
    cells = colnames(obj)
  )

  cat(glue(
    "[{smp}] FeatureMatrix completed.\n",
    "  Condition: {cond}\n",
    "  Matrix dim: {nrow(mat)} x {ncol(mat)}\n"
  ))

  cat(glue("[{smp}] Inserting unified peaks assay into SeuratObject...\n"))
  obj[["peaks"]] <- CreateChromatinAssay(
    counts = mat,
    fragments = Fragments(obj)[[1]],
    annotation = annotation
  )

  cat(glue(
    "[{smp}] Unified peaks assay inserted.\n",
    "  Assays: {paste(Assays(obj), collapse = ', ')}\n",
    "  peaks features: {nrow(obj[['peaks']])}\n"
  ))

  out_rds <- file.path(seurat_dir, paste0(smp, "_qc_noMT_macs3.rds"))
  cat(glue("[{smp}] Saving SeuratObject with unified peaks assay...\n"))
  saveRDS(obj, out_rds)
  cat(glue("[{smp}] Saved in: {out_rds}\n\n"))
}

for (i in seq_len(nrow(samples))) {
  smp <- samples$sample[i]
  cond <- samples$condition[i]

  obj_file <- file.path(seurat_dir, paste0(smp, "_qc_noMT_macs3.rds"))
  cat(glue("[{smp}] Loading SeuratObject with unified peaks assay from: {obj_file}\n"))
  obj <- readRDS(obj_file)

  export_dir <- file.path(macs_out_root, paste0(smp, "_peaks_matrix"))
  dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

  cat(glue("[{smp}] Exporting peaks assay in Python-readable format...\n"))

  peak_counts <- LayerData(obj, assay = "peaks", layer = "counts")

  # Save sparse matrix as Matrix Market
  mtx_file <- file.path(export_dir, "matrix.mtx")
  Matrix::writeMM(obj = peak_counts, file = mtx_file)
  system2("gzip", args = c("-f", mtx_file))

  # Save barcodes
  barcodes_file <- file.path(export_dir, "barcodes.tsv.gz")
  writeLines(colnames(peak_counts), gzfile(barcodes_file))

  # Save features as chr:start-end
  feature_names <- rownames(peak_counts)
  features_file <- file.path(export_dir, "features.tsv.gz")
  writeLines(feature_names, gzfile(features_file))

  # Save peaks as BED-like table
  peak_gr <- StringToGRanges(feature_names, sep = c(":", "-"))
  peaks_bed <- data.frame(
    chrom = as.character(seqnames(peak_gr)),
    start = start(peak_gr) - 1,
    end = end(peak_gr),
    peak_id = feature_names,
    stringsAsFactors = FALSE
  )
  peaks_bed_file <- file.path(export_dir, "peaks.bed.gz")
  write.table(
    peaks_bed,
    file = gzfile(peaks_bed_file),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # Save fragments path used in this object
  frag_obj <- Fragments(obj)[[1]]
  frag_path <- frag_obj@path
  frag_info_file <- file.path(export_dir, "fragments_path.txt")
  writeLines(frag_path, frag_info_file)

  cat(glue(
    "[{smp}] Python-readable peaks matrix exported.\n",
    "  Condition: {cond}\n",
    "  Matrix: {file.path(export_dir, 'matrix.mtx.gz')}\n",
    "  Barcodes: {barcodes_file}\n",
    "  Features: {features_file}\n",
    "  Peaks BED: {peaks_bed_file}\n",
    "  Fragments path record: {frag_info_file}\n\n"
  ))
}

cat(glue(
  "[DONE] MACS3 peak calling, union peak generation, unified peaks assay insertion, and Python-readable peaks matrix export completed for all samples.\n"
))
