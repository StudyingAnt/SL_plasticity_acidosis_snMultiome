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
MEM=200G TIME=12:00:00 bash slurm/submit.sh 00-2_seurat_preprocessing/05_transformations.R
```

Logs are written to `LOG_DIR` as `snmult_<step>.<jobid>.out/.err`. Each log records the repository commit that was run.
