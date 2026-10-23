# Single-nuclei multi-omics analysis of sphingolipid metabolic plasticity under acidosis

## Introduction

This repository contains the scripts for the single-nucleus multi-omics analysis (10x Genomics Multiome, joint snRNA-seq + snATAC-seq) in the paper below: Figure 6 and the corresponding Methods. MCF7 cells were profiled under four conditions:

| Sample | Condition | Description |
|---|---|---|
| ct74 | Control | non-adapted cells, pH 7.4 |
| ct65 | Acute | non-adapted cells, pH 6.5 for 3 days |
| aa65 | Chronic | acid-adapted cells, pH 6.5 for > 3 months |
| aa74 | Return | acid-adapted cells returned to pH 7.4 for 3 days |

## Paper

Chalar R., Khatri N., Obeid J., Song J.-H., *et al.* Metabolic plasticity enables cancer adaptation to acidosis in tumor ecosystems. *bioRxiv* (2026). [doi:10.64898/2026.01.14.699600](https://doi.org/10.64898/2026.01.14.699600). See [CITATION.cff](CITATION.cff).

## Start here

1. **Data:** raw FASTQ and processed Cell Ranger ARC files are in GEO ([DATA_AVAILABILITY.md](DATA_AVAILABILITY.md)).
2. **Environment:** see [env/README.md](env/README.md).
3. **Configuration:** copy `config/config.example.R` to `config/config.R` and `config/config.example.sh` to `config/config.sh`, then edit the paths. Keep results outside the repository.
4. **Run:** run the steps below in order, from the repository root. Each folder's README lists inputs, outputs, parameters and expected results.

## Pipeline

| Step | Folder | What it does | Figure |
|---|---|---|---|
| 0-1 | [00-1_cellranger_arc](00-1_cellranger_arc) | FASTQ → Cell Ranger ARC 2.0.2 (GRCh38) | – |
| 0-2 | [00-2_seurat_preprocessing](00-2_seurat_preprocessing) | Seurat/Signac objects, QC, MACS3 peaks, SCT/LSI/WNN, ACDC clustering, chromVAR | – |
| 1 | [01_metacells](01_metacells) | hdWGCNA metacells (k = 25), metacell CPM matrix | – |
| 2 | [02_aracne3](02_aracne3) | Human TF list, TF regulatory network with ARACNe3 | – |
| 3 | [03_viper](03_viper) | Regulon, VIPER activity, metacell UMAPs, metacell chromVAR | 6a |
| 4 | [04_sphingolipid_axes](04_sphingolipid_axes) | Regulatory pressure on four ceramide clearance axes, dominant axis, metabolic diversity | 6b–e |
| 5 | [05_axis_regulators](05_axis_regulators) | msVIPER candidate regulators per axis with network linkage and motif support | 6f–i |
| 6 | [06_extended_data_atac](06_extended_data_atac) | RNA/ATAC/WNN UMAPs by condition, chromVAR deviations of the candidate regulators | Extended Data |

Example:

```bash
bash 00-1_cellranger_arc/run_cellranger_arc.sh
Rscript 00-2_seurat_preprocessing/01_make_objects.R
...
Rscript 05_axis_regulators/03_plot_regulator_dotplots.R
```

On a SLURM cluster, the whole pipeline can be submitted as one dependency chain with `bash slurm/run_all.sh` (see [slurm/README.md](slurm/README.md)).

## Reproducibility

The full pipeline was rerun from the Cell Ranger ARC outputs in a separate environment. QC, clustering and metacells were identical to the original analysis, the ARACNe3 network contained all original edges, and every Figure 6 panel was reproduced with the same candidate regulators. Details, and the differences from the original notebooks, are in [NOTES.md](NOTES.md).

## Repository layout

```
00-1_cellranger_arc/      Cell Ranger ARC (samples.csv, libraries/, run script, author_hpc/)
00-2_seurat_preprocessing/
01_metacells/  02_aracne3/  03_viper/  04_sphingolipid_axes/  05_axis_regulators/  06_extended_data_atac/
R/common.R                shared settings (samples, colors, axis gene sets, output folders)
config/                   config templates and QC thresholds
slurm/                    SLURM submission (submit.sh, run_all.sh)
env/                      software versions and environment files
```

## License

[MIT](LICENSE)
