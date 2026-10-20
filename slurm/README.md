# Running on a SLURM cluster

`submit.sh` submits one pipeline step as a SLURM job. `job_step.sh` is the job body: it activates the R environment, runs the step from the repository root, and prints `sessionInfo()` at the end.

Setup (once):

1. Copy `config/config.example.R` to `config/config.R`, and `config/config.example.sh` to `config/config.sh`, and fill in the paths.
2. In `config/config.sh`, set `CONDA_MODULE`, `CONDA_ENV` and `R_LIBRARY` (the R environment), plus `SLURM_PARTITION`, `LOG_DIR` and optionally `SLURM_EXTRA_ARGS`.

Submit a step:

```bash
bash slurm/submit.sh 00-2_seurat_preprocessing/01_make_objects.R
```

Each step has default resources (CPUs, memory, time) in `submit.sh`. To override them:

```bash
MEM=200G TIME=12:00:00 PARTITION=long-96core-shared bash slurm/submit.sh 00-2_seurat_preprocessing/05_transformations.R
```

Logs are written to `LOG_DIR` as `snmult_<step>.<jobid>.out/.err`. Each log records the repository commit that was run, the config files and the versions of the main R packages.

To run with other config files (for example a test run with its own `RESULTS_DIR`), set `SNMULT_CONFIG_R` and `SNMULT_CONFIG_SH`. `JOB_TAG` is added to the job and log names. Any `config/config*.R` and `config/config*.sh` file is ignored by git.

```bash
SNMULT_CONFIG_R=config/config.test.R SNMULT_CONFIG_SH=config/config.test.sh JOB_TAG=test \
  bash slurm/submit.sh 01_metacells/01_build_metacells.R
```

To chain steps, pass the job id of the previous step with `AFTER`; the next job starts only if that job finished successfully. `submit.sh` prints the new job id on its last line:

```bash
j1=$(bash slurm/submit.sh 03_viper/03_metacell_embeddings.R | tail -1)
AFTER=$j1 bash slurm/submit.sh 03_viper/04_attach_chromvar.R
```
