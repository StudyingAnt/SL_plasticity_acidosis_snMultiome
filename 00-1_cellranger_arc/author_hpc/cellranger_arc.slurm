#!/bin/bash
#SBATCH --job-name=cellranger_arc
#SBATCH --output=/gpfs/projects/DamaghiGroup/projects/mcf7_aa_cacr_lipid/logs/cellranger_arc.%A_%a.out
#SBATCH --error=/gpfs/projects/DamaghiGroup/projects/mcf7_aa_cacr_lipid/logs/cellranger_arc.%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=24:00:00
#SBATCH --partition=long-40core-shared
#SBATCH --mail-user=<your-email>
#SBATCH --mail-type=END,FAIL

set -euo pipefail

PROJECT=${PROJECT:-/gpfs/projects/DamaghiGroup/projects/mcf7_aa_cacr_lipid}
REF=/gpfs/projects/DamaghiGroup/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
ARC=/gpfs/projects/DamaghiGroup/software/cellranger-arc-2.0.2/cellranger-arc
SCRATCH_ROOT=/gpfs/scratch/joonsong/cellranger_arc

: "${SLURM_ARRAY_TASK_ID:?Submit as an array job}"
: "${SLURM_CPUS_PER_TASK:?Missing CPU allocation}"
: "${SLURM_MEM_PER_NODE:?Missing memory allocation}"

SAMPLE=$(python3 - "$PROJECT/config/samples.csv" "$SLURM_ARRAY_TASK_ID" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as handle:
    rows = list(csv.DictReader(handle))
idx = int(sys.argv[2])
if not 0 <= idx < len(rows):
    raise SystemExit(f'Invalid array index {idx}')
print(rows[idx]['sample_id'])
PY
)
LIBRARIES="$PROJECT/config/libraries/$SAMPLE.csv"
DEST="$PROJECT/data/cellranger_arc/$SAMPLE"
WORK="$SCRATCH_ROOT/${SAMPLE}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
META="$PROJECT/runs/cellranger_arc/${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "$META" "$PROJECT/data/cellranger_arc" "$SCRATCH_ROOT"

if [[ -e "$DEST" ]]; then
    echo "Destination already exists; refusing overwrite: $DEST" >&2
    exit 2
fi
for target in "$ARC" "$REF/reference.json" "$LIBRARIES"; do
    [[ -e "$target" ]] || { echo "Missing input $target" >&2; exit 2; }
done
[[ ! -e "$WORK" ]] || { echo "Scratch directory already exists: $WORK" >&2; exit 2; }
mkdir "$WORK"

# Keep the raw FASTQ files in their shared original storage location.
# Calculate with a modest memory margin below the SLURM allocation.
ALLOC_MB="$SLURM_MEM_PER_NODE"
LOCALMEM=$((ALLOC_MB / 1024 - 8))
((LOCALMEM > 0)) || { echo 'Insufficient requested memory' >&2; exit 2; }
{
    date -Is
    printf 'sample=%s\njob_id=%s\nnode=%s\ncpus=%s\nslurm_memory_mb=%s\narc_localmem_gb=%s\nwork=%s\n' \
        "$SAMPLE" "$SLURM_JOB_ID" "$(hostname)" "$SLURM_CPUS_PER_TASK" "$ALLOC_MB" "$LOCALMEM" "$WORK"
    "$ARC" --version
} > "$META/run_info.txt" 2>&1
cp "$LIBRARIES" "$META/libraries.csv"
cd "$WORK"
set +e
"$ARC" count --id="$SAMPLE" --reference="$REF" --libraries="$LIBRARIES" \
    --localcores="$SLURM_CPUS_PER_TASK" --localmem="$LOCALMEM" \
    > "$META/cellranger_stdout.log" 2> "$META/cellranger_stderr.log"
status=$?
set -e
printf '%s\n' "$status" > "$META/exit_code.txt"
if ((status != 0)); then
    echo "Cell Ranger ARC failed ($status). Scratch preserved: $WORK" >&2
    exit "$status"
fi
[[ -d "$WORK/$SAMPLE/outs" ]] || { echo "Missing outs directory; scratch preserved: $WORK" >&2; exit 3; }

# Stage results to a hidden directory before exposing the final sample path.
STAGE="$PROJECT/data/cellranger_arc/.${SAMPLE}_${SLURM_JOB_ID}.staging"
[[ ! -e "$STAGE" ]] || { echo "Staging already exists: $STAGE" >&2; exit 3; }
mkdir "$STAGE"
if ! rsync -a -- "$WORK/$SAMPLE/outs/" "$STAGE/outs/"; then
    echo "Transfer failed; stage and scratch preserved" >&2
    exit 4
fi
if ! rsync -a --checksum --dry-run --itemize-changes -- "$WORK/$SAMPLE/outs/" "$STAGE/outs/" > "$META/verification.txt"; then
    echo "Verification failed; stage and scratch preserved" >&2
    exit 4
fi
if [[ -s "$META/verification.txt" ]]; then
    echo "Verification found differences; stage and scratch preserved" >&2
    exit 4
fi
[[ ! -e "$DEST" ]] || { echo "Destination appeared during transfer; scratch preserved" >&2; exit 4; }
mv -T "$STAGE" "$DEST"
printf 'success\t%s\n' "$(date -Is)" > "$META/status.txt"
# Delete only this task's owned scratch directory after verified transfer.
rm -rf -- "$WORK"
echo "Completed: $SAMPLE -> $DEST"
