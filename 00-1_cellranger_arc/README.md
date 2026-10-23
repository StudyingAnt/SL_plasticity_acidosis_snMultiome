# Cell Ranger ARC

This step aligns and counts the raw 10x Multiome reads (gene expression + ATAC) for the four samples with Cell Ranger ARC. It produces the per-sample matrices and fragment files that [00-2_seurat_preprocessing](../00-2_seurat_preprocessing) reads.

| Sample | Condition | GEX FASTQ | ATAC FASTQ |
|---|---|---|---|
| ct74 | Control: non-adapted MCF7 at pH 7.4 | SC60_1G | SC60_1A |
| ct65 | Acute: non-adapted MCF7 at pH 6.5 for 3 days | SC60_2G | SC60_2A |
| aa65 | Chronic: acid-adapted MCF7 at pH 6.5 (> 3 months) | SC60_4G | SC60_4A |
| aa74 | Return: acid-adapted MCF7 returned to pH 7.4 for 3 days | SC60_3G | SC60_3A |

The full sample sheet, including the FASTQ name prefixes, is [samples.csv](samples.csv).

## Requirements

- **Cell Ranger ARC 2.0.2**, downloaded from the 10x Genomics website.
- **Reference `refdata-cellranger-arc-GRCh38-2020-A-2.0.0`** (10x Genomics GRCh38 reference), downloaded with:

```bash
curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
tar xzf refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz
```

- **Raw FASTQ files** from GEO ([DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md)). Each library has I1, R1 and R2 (GEX) or I1, R1, R2 and R3 (ATAC) files. Cell Ranger ARC requires Illumina file names (`<prefix>_S<n>_L<lane>_<read>_001.fastq.gz`). The original names are listed in the GEO records and in [samples.csv](samples.csv).
- About 64 GB of memory and 16 cores per sample. Each sample took 8–9 hours with these resources.

## Run

1. Put the FASTQ files in one folder per library:

```
<FASTQ_DIR>/scRNA/SC60_1G/SC60_1G_CKDL230025675-1A_225YMKLT3_S3_L006_R1_001.fastq.gz ...
<FASTQ_DIR>/scATAC/SC60_1A/SC60_1A_CKDL230025679-1A_225VWWLT3_S4_L008_R1_001.fastq.gz ...
```

2. Set `CELLRANGER_ARC`, `ARC_REFERENCE`, `FASTQ_DIR`, `ARC_OUT_DIR`, `ARC_CORES` and `ARC_MEM_GB` in `config/config.sh`.

3. Run all four samples, or one sample at a time:

```bash
bash 00-1_cellranger_arc/run_cellranger_arc.sh
bash 00-1_cellranger_arc/run_cellranger_arc.sh ct74
```

For each sample, the script writes a libraries CSV to `<ARC_OUT_DIR>/libraries/` and runs:

```bash
cellranger-arc count --id=<sample> --reference=<ARC_REFERENCE> --libraries=<sample>.csv \
  --localcores=16 --localmem=64
```

The libraries CSVs used for the paper are in [libraries/](libraries) (with `<FASTQ_DIR>` in place of the local path).

## Output

`<ARC_OUT_DIR>/<sample>/outs/`. Set `DATASET_DIR` in `config/config.R` to `ARC_OUT_DIR`. The next step uses two of the output files:

| File | Used for |
|---|---|
| `filtered_feature_bc_matrix.h5` | gene expression and ATAC peak counts of cell-associated barcodes |
| `atac_fragments.tsv.gz` (+ `.tbi`) | ATAC fragments (QC metrics, MACS3 peak calling, peak counting) |

## Expected results

From `outs/summary.csv`. The full metrics are in [author_hpc/run_records/cellranger_arc_summary.csv](author_hpc/run_records/cellranger_arc_summary.csv).

| Sample | Estimated cells | GEX read pairs | GEX median genes / cell | ATAC read pairs | ATAC median fragments / cell | ATAC TSS enrichment | ATAC peaks |
|---|---|---|---|---|---|---|---|
| ct74 | 10,371 | 293,796,227 | 3,899 | 306,005,271 | 19,425 | 10.66 | 203,621 |
| ct65 | 13,810 | 249,077,125 | 3,096 | 298,098,905 | 14,041 | 10.82 | 217,992 |
| aa65 | 10,464 | 265,707,135 | 3,855 | 298,180,391 | 18,826 | 11.12 | 200,804 |
| aa74 | 10,225 | 299,868,478 | 4,030 | 337,076,440 | 19,768 | 11.21 | 197,769 |

After QC in 00-2, these give 9,318 / 12,774 / 9,821 / 9,442 nuclei (ct74 / ct65 / aa65 / aa74). These are the same numbers as in the analysis for the paper.

## Author execution

The authors ran this step on the SeaWulf cluster at Stony Brook University (SLURM, Cell Ranger ARC 2.0.2, Martian 4.0.5, STAR 2.7.2a). The scripts and run records are in [author_hpc/](author_hpc):

| File | Description |
|---|---|
| `setup_cellranger_arc.py` | Links the raw FASTQ folders into the project without modifying them, and writes one libraries CSV per sample |
| `submit_cellranger_arc.sh` | Submits the four samples as one SLURM array (two at a time, 16 cores, 72 GB) |
| `cellranger_arc.slurm` | Runs `cellranger-arc count` in scratch, copies `outs/` back with checksum verification, and records the run |
| `samples.csv` | Sample sheet used on the cluster |
| `run_records/` | Job, node, resources, start/end time, exit code and transfer check for each sample (`cellranger_arc_runs.tsv`); tool versions; `summary.csv` metrics |

Paths in these files are specific to SeaWulf; adapt them to your site.
