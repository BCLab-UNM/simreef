#!/bin/bash

# Default cluster name
CLUSTER=""

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --cluster)
            CLUSTER="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            echo "Usage: $0 --cluster <clustername>"
            exit 1
            ;;
    esac
done

# Require cluster name
if [ -z "$CLUSTER" ]; then
    echo "Error: --cluster argument is required"
    echo "Usage: $0 --cluster <clustername>"
    exit 1
fi

cd slurm || { echo "Failed to cd into slurm directory"; exit 1; }

SLURM_SCRIPT="debug_${CLUSTER}.slurm"

# Submit the job and capture Job ID
submit_output=$(sbatch "$SLURM_SCRIPT")
if [ $? -ne 0 ]; then
  echo "sbatch failed: $submit_output"
  exit 1
fi

jobid=$(echo "$submit_output" | awk '{print $4}')
echo "Submitted job $jobid, waiting for matching output file(s)..."

# Wait for any file containing the job ID
while true; do
  matches=$(grep -l "$jobid" * 2>/dev/null)
  if [ -n "$matches" ]; then
    break
  fi
  sleep 1
done

# Display output files
echo $PWD
echo -e "\nOutput files containing job ID $jobid:"
echo "$matches"

echo -e "\n---- SLURM OUTPUT BEGIN ----"
tail -f "$matches"
