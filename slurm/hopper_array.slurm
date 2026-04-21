#!/bin/bash
#SBATCH --account=2016345
#SBATCH --partition=biocomp,general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --job-name=simreef_params
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=64GB
#SBATCH --array=1-1

set -euo pipefail

# Setup environment
export UPCXX_SHARED_HEAP_SIZE='2 GB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv

# GASNet settings
export GASNET_MAX_SEGSIZE='2 GB'
export GASNET_PHYSMEM_MAX='2 GB'
export GASNET_IBV_PHYSMEM_MAX='2 GB'
export GASNET_USE_ODP=0
export GASNET_BACKTRACE=1

# Load modules
module load gcc/12.1.0-crtl
module load upcxx/2020.10.0-6eh2
module load opencv/4.10.0-uddb
module load python

SIMFORAGER_PATH="$HOME/simreef/simforager/install/bin/simforager"
CONFIG_TEMPLATE_FILE="$HOME/simreef/simforager/configs/simforager_template.config"
PARAMETER_TABLE="$HOME/simreef/simforager/configs/parameter_table.csv"

if [ ! -e "$SIMFORAGER_PATH" ]; then
    echo "The executable $SIMFORAGER_PATH does not exist. Exiting."
    exit 1
fi

if [ ! -e "$CONFIG_TEMPLATE_FILE" ]; then
    echo "The configuration template $CONFIG_TEMPLATE_FILE does not exist. Exiting."
    exit 1
fi

if [ ! -e "$PARAMETER_TABLE" ]; then
    echo "The parameter table $PARAMETER_TABLE does not exist. Exiting."
    exit 1
fi

NUM_ROWS=$(($(wc -l < "$PARAMETER_TABLE") - 1))
if [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "$NUM_ROWS" ]; then
    echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is out of range. Table has $NUM_ROWS data rows."
    exit 1
fi

ROW_ID=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")

OUTPUT_DIR="/carc/scratch/users/mfricke/simreef_outputs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}/${ROW_ID}"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

CONFIG_FILE="$OUTPUT_DIR/simforager.config"
cp "$CONFIG_TEMPLATE_FILE" "$CONFIG_FILE"

# Extract the selected CSV row as key/value pairs.
python3 - "$PARAMETER_TABLE" "$SLURM_ARRAY_TASK_ID" > "$OUTPUT_DIR/selected_params.tsv" <<'PY'
import csv
import sys

table_path = sys.argv[1]
row_index = int(sys.argv[2])

with open(table_path, newline='') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

row = rows[row_index - 1]

for key, value in row.items():
    print(f"{key}\t{value}")
PY

# Replace placeholders in the config template using sed.
while IFS=$'\t' read -r key value; do
    placeholder=$(echo "$key" | tr '[:lower:]-' '[:upper:]_')
    escaped_value=$(printf '%s\n' "$value" | sed 's/[\/&]/\\&/g')
    sed -i "s|<${placeholder}>|${escaped_value}|g" "$CONFIG_FILE"
done < "$OUTPUT_DIR/selected_params.tsv"

echo "Selected parameter row: $SLURM_ARRAY_TASK_ID"
echo "Output directory: $OUTPUT_DIR"
echo "Final config file: $CONFIG_FILE"

CMD="srun --mpi=pmi2 $SIMFORAGER_PATH --config $CONFIG_FILE"
echo "Running: $CMD"

srun --mpi=pmi2 "$SIMFORAGER_PATH" --config "$CONFIG_FILE"
