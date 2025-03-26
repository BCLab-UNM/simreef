#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --job-name=simreef
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=24GB
#SBATCH --output=debug_%j.out
#SBATCH --error=debug_%j.out

# Setup environment
export UPCXX_SHARED_HEAP_SIZE='2 GB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv  # Try 'mpi' if OOM persists

# GASNet settings
export GASNET_MAX_SEGSIZE='2 GB'
export GASNET_PHYSMEM_MAX='2 GB'
export GASNET_IBV_PHYSMEM_MAX='2 GB'
export GASNET_USE_ODP=0  # Disable On-Demand Paging if OOM occurs

# Print backtrace on crash
export GASNET_BACKTRACE=1

# Load modules
module load gcc/12.1.0-crtl
module load upcxx/2020.10.0-6eh2

SIMFORAGER_PATH=$HOME/simreef/simforager/install/bin/simforager
CONFIG_FILE=$HOME/simreef/simforager/configs/simreef_visual.config

# Check that binary and config both exist
if [ ! -e "$SIMFORAGER_PATH" ]; then
    echo "The executable $SIMFORAGER_PATH does not exist. Exiting."
    exit 1
fi

if [ ! -e "$CONFIG_FILE" ]; then
    echo "The configuration file $CONFIG_FILE does not exist. Exiting."
    exit 1
fi

# Create an output path based on the job name
OUTPUT_DIR=$HOME/simreef/outputs/${SLURM_JOB_NAME}_${SLURM_JOBID}_$(date +'%Y-%m-%d_%H-%M-%S')/
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo Running: srun --mpi=pmi2 $SIMFORAGER_PATH --config $CONFIG_FILE in $OUTPUT_DIR

srun --mpi=pmi2 $SIMFORAGER_PATH --config $CONFIG_FILE
