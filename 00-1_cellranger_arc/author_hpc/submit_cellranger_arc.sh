#!/usr/bin/env bash
set -euo pipefail
PROJECT=${PROJECT:-/gpfs/projects/DamaghiGroup/projects/mcf7_aa_cacr_lipid}
PARTITION=${PARTITION:-long-40core-shared}
CONCURRENCY=${CONCURRENCY:-2}
CPUS=${CPUS:-16}
MEM=${MEM:-72G}
TIME=${TIME:-24:00:00}

if (($#)); then
    echo 'Usage: PARTITION=... CONCURRENCY=2 bash code/submit_cellranger_arc.sh' >&2
    exit 2
fi
python3 "$PROJECT/code/python/setup_cellranger_arc.py" --project "$PROJECT"
mkdir -p "$PROJECT/logs" "$PROJECT/runs/cellranger_arc" "$PROJECT/data/cellranger_arc"
count=$(python3 - "$PROJECT/config/samples.csv" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    print(sum(1 for _ in csv.DictReader(f)))
PY
)
(( count > 0 )) || { echo 'No samples' >&2; exit 2; }
[[ "$CONCURRENCY" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid CONCURRENCY' >&2; exit 2; }
# Run a single SLURM array with bounded concurrent jobs.
sbatch --export="ALL,PROJECT=$PROJECT" --partition="$PARTITION" \
    --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME" \
    --array="0-$((count-1))%${CONCURRENCY}" \
    "$PROJECT/slurm/cellranger_arc.slurm"
