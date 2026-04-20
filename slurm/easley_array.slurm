#!/bin/bash
#SBATCH --account=2016345
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --job-name=simreef_easley
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=64G
#SBATCH --array=1-64

set -euo pipefail

# Make sure the storage group membership is active
newgrp melaniem2016345

# UPC++ / GASNet environment
export UPCXX_SHARED_HEAP_SIZE='2 GB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=debug
export UPCXX_NETWORK=ibv

export GASNET_MAX_SEGSIZE='2 GB'
export GASNET_PHYSMEM_MAX='2 GB'
export GASNET_IBV_PHYSMEM_MAX='2 GB'
export GASNET_USE_ODP=0
export GASNET_BACKTRACE=1

# Load Easley build-compatible modules
module purge
module load boost/1.88.0-65dn
module load openmpi/4.1.7-e7k3
module load upcxx/2023.9.0-27k4
module load cmake/3.31.6-auvg
module load opencv/4.10.0-ki66
module load libdeflate/1.14-2pby
module load python

SIMFORAGER_PATH="$HOME/simreef/simforager/install_easley/bin/simforager"
CONFIG_TEMPLATE_FILE="$HOME/simreef/simforager/configs/simforager_template.config"
PARAMETER_TABLE="$HOME/simreef/simforager/configs/parameter_table.csv"

export LD_LIBRARY_PATH="${OPENCV_LIB64:-}:${OPENCV_LIB:-}:${CARC_LIBRARY_PATH:-}:${LD_LIBRARY_PATH:-}"fi

if [ ! -x "$SIMFORAGER_PATH" ]; then
    echo "Executable not found or not executable: $SIMFORAGER_PATH"
    exit 1
fi

if [ ! -f "$CONFIG_TEMPLATE_FILE" ]; then
    echo "Configuration template not found: $CONFIG_TEMPLATE_FILE"
    exit 1
fi

if [ ! -f "$PARAMETER_TABLE" ]; then
    echo "Parameter table not found: $PARAMETER_TABLE"
    exit 1
fi

NUM_ROWS=$(($(wc -l < "$PARAMETER_TABLE") - 1))
if [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "$NUM_ROWS" ]; then
    echo "SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID is out of range. Table has $NUM_ROWS data rows."
    exit 1
fi

ROW_ID=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")

OUTPUT_BASE="/carc/scratch/users/mfricke/simreef_outputs"
OUTPUT_DIR="${OUTPUT_BASE}/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}/${ROW_ID}"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

CONFIG_FILE="$OUTPUT_DIR/simforager.config"
cp "$CONFIG_TEMPLATE_FILE" "$CONFIG_FILE"

# Extract the selected CSV row as key/value pairs
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

# Replace placeholders in config template
while IFS=$'\t' read -r key value; do
    placeholder=$(echo "$key" | tr '[:lower:]-' '[:upper:]_')
    escaped_value=$(printf '%s\n' "$value" | sed 's/[\/&]/\\&/g')
    sed -i "s|<${placeholder}>|${escaped_value}|g" "$CONFIG_FILE"
done < "$OUTPUT_DIR/selected_params.tsv"

echo "=================================================="
echo "Job name:          $SLURM_JOB_NAME"
echo "Job ID:            $SLURM_JOB_ID"
echo "Array job ID:      $SLURM_ARRAY_JOB_ID"
echo "Array task ID:     $SLURM_ARRAY_TASK_ID"
echo "Selected row:      $ROW_ID"
echo "Output directory:  $OUTPUT_DIR"
echo "Executable:        $SIMFORAGER_PATH"
echo "Config file:       $CONFIG_FILE"
echo "=================================================="

echo "Resolved parameters:"
cat "$OUTPUT_DIR/selected_params.tsv"
echo "=================================================="

# Avoid UCX PML on Easley. UCX is crashing during MPI_Init.
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp
export OMPI_MCA_osc=pt2pt

CMD=(srun --mpi=pmi2 "$SIMFORAGER_PATH" --config "$CONFIG_FILE")
echo "Running: ${CMD[*]}"

"${CMD[@]}"
