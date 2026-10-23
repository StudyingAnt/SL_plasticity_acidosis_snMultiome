# Software environment

The whole pipeline was rerun and verified in the environment below (SLURM cluster, Rocky Linux 9; see [NOTES.md](../NOTES.md)).

## External tools

| Tool | Version | Step | Notes |
|---|---|---|---|
| Cell Ranger ARC | 2.0.2 | 00-1 | 10x Genomics; reference `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` |
| MACS3 | 3.0.3 (Python 3.10.19) | 00-2 step 03 | [macs3.yml](macs3.yml) |
| ARACNe3 | califano-lab/ARACNe3 (commit `3d8791a` for the original analysis) | 02 | Built from source (https://github.com/califano-lab/ARACNe3); run without a cached null model (see NOTES.md) |

## R

R 4.5.2 from conda ([r452.yml](r452.yml)); R packages from Bioconductor 3.22 / CRAN, managed with renv ([renv.lock](renv.lock)).

| Package | Version |
|---|---|
| Seurat / SeuratObject | 5.4.0 / 5.3.0 |
| Signac | 1.16.0 |
| sctransform | 0.4.3 (glmGamPoi **not** installed, so SCTransform v2 uses `method = "nb_offset"` as in the original analysis) |
| leidenbase | 0.1.37 |
| chromVAR / motifmatchr | 1.32.0 / 1.32.0 |
| JASPAR2024 / TFBSTools | 0.99.7 / 1.48.0 |
| EnsDb.Hsapiens.v86 / BSgenome.Hsapiens.UCSC.hg38 | 2.99.0 / 1.4.5 |
| hdWGCNA | 0.4.12 |
| viper / genefilter | 1.44.0 / 1.92.0 |
| data.table / ggplot2 / Matrix | 1.18.2.1 / 4.0.2 / 1.7-4 |

## Setup

Create the environments outside the repository (e.g. in `/path/to/envs`):

```bash
conda env create -p /path/to/envs/r452 -f env/r452.yml
conda env create -p /path/to/envs/macs3 -f env/macs3.yml
mkdir -p /path/to/envs/r452_renv && cp env/renv.lock /path/to/envs/r452_renv/
conda activate /path/to/envs/r452
cd /path/to/envs/r452_renv
Rscript -e 'renv::init(bare = TRUE, restart = FALSE); renv::restore(lockfile = "renv.lock", prompt = FALSE)'
Rscript -e 'source("renv/activate.R"); cat(.libPaths()[1], "\n")'   # -> R_LIBRARY in config/config.sh
```

`renv.lock` records 304 packages (R 4.5.2, Bioconductor 3.22), including the hdWGCNA GitHub commit.

If an R package fails to compile against conda headers, pointing `CPATH` at the conda include directories helps, e.g. `export CPATH=$CONDA_PREFIX/include/freetype2:$CONDA_PREFIX/include`.

## Original analysis

Steps 00-2 were originally run on the same cluster with R 4.5.2 and Bioconductor 3.22. Steps 01–05 were originally run on a Windows workstation with R 4.3.2 and hdWGCNA 0.4.10, and ARACNe3 was run under WSL. The rerun in the environment above reproduced the original results (see [NOTES.md](../NOTES.md)).
