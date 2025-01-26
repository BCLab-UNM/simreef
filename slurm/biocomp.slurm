#!/bin/bash
#SBATCH --partition=biocomp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00-30:00
#SBATCH --job-name=simreef
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL

# Setup environment
export GASNET_MAX_SEGSIZE=128MB/P
export UPCXX_SHARED_HEAP_SIZE='128 MB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv

# Show warnings
GASNET_ODP_VERBOSE=1

# Print backtrace on crash
GASNET_BACKTRACE=1

# For Biocomp we can assume 124 GB of RAM
GASNET_PHYSMEM_MAX='124 GB'

# Load modules

module load gcc/12.1.0-crtl
module load upcxx/2020.10.0-6eh2

SIMFORAGER_PATH=$HOME/simreef/simforager/install/bin/simforager
CONFIG_FILE=$HOME/simreef/simforager/configs/simforager_test.config

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

