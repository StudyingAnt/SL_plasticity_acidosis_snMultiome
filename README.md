# Single-nuclei multi-omics analysis of sphingolipid metabolic plasticity under acidosis

## Introduction

This repository contains the scripts used for the single-nuclei multi-omics (10x Multiome, snRNA + snATAC) analysis in the paper below (Figure 6 and the corresponding Methods sections). MCF7 cells were profiled under four conditions:

| Sample | Condition | Description |
|---|---|---|
| ct74 | Control | non-adapted cells, pH 7.4 |
| ct65 | Acute | non-adapted cells, pH 6.5 for 3 days |
| aa65 | Chronic | acid-adapted cells, pH 6.5 for > 3 months |
| aa74 | Return | acid-adapted cells returned to pH 7.4 for 3 days |

## Paper

Chalar R., Khatri N., Obeid J., Song J.-H., *et al.* Metabolic plasticity enables cancer adaptation to acidosis in tumor ecosystems. *(under revision)*

## Data availability

Raw and processed sequencing data: *(accession to be added)*

## Pipeline

Each folder has its own README with inputs, outputs and parameters. Run the steps in order.

| Step | Folder | What it does | Figure |
|---|---|---|---|
| 0-1 | [00-1_cellranger_arc](00-1_cellranger_arc) | FASTQ → Cell Ranger ARC (GRCh38) | – |
| 0-2 | [00-2_seurat_preprocessing](00-2_seurat_preprocessing) | Seurat/Signac objects, QC, MACS3 peaks, SCT/LSI/WNN, ACDC clustering, chromVAR | – |
| 1 | [01_metacells](01_metacells) | hdWGCNA metacells (k = 25), metacell CPM matrix | – |
| 2 | [02_aracne3](02_aracne3) | TF regulatory network with ARACNe3 | – |
| 3 | [03_viper](03_viper) | Regulon, VIPER activity, metacell UMAPs, metacell chromVAR | 6a |
| 4 | [04_sphingolipid_axes](04_sphingolipid_axes) | Regulatory pressure on four ceramide clearance axes, dominant axis, metabolic diversity | 6c–f |
| 5 | [05_axis_regulators](05_axis_regulators) | msVIPER candidate regulators per axis with network linkage and motif support | 6g–j |

## How to run

1. Install the software listed in [env/README.md](env/README.md).
2. Copy `config/config.example.R` to `config/config.R` and `config/config.example.sh` to `config/config.sh`, then edit the paths.
3. Run every script from the repository root, for example:

```bash
Rscript 00-2_seurat_preprocessing/01_make_objects.R
Rscript 00-2_seurat_preprocessing/06_clustering_acdc.R -i /path/to/seurat_objects/all_qc_noMT_macs3_trf.rds
Rscript 01_metacells/01_build_metacells.R
bash 02_aracne3/02_run_aracne3.sh
```

On a SLURM cluster, steps can be submitted with `bash slurm/submit.sh <step script>` (see [slurm/README.md](slurm/README.md)).

Step 0-2 was run on an HPC cluster (Linux). Steps 1–5 were run on a workstation; ARACNe3 was run under Linux (WSL).

## Notes

The scripts were reorganized from the analysis notebooks used for the paper. Paths were replaced by `config/config.R`, and exploratory code that did not contribute to the paper was removed. Analysis parameters were not changed. Remaining differences are listed in [NOTES.md](NOTES.md).

## License

[MIT](LICENSE)
