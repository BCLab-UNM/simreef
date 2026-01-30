#!/bin/bash
#SBATCH --account 2016345
#SBATCH --partition=condo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --job-name=simreef_array
#SBATCH --mail-user=htasnim30@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=24GB
#SBATCH --array 1-30

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
module load opencv/4.10.0-uddb

SIMFORAGER_PATH=$HOME/simreef/simforager/install/bin/simforager
CONFIG_TEMPLATE_FILE=$HOME/simreef/simforager/configs/simforager_template_social.config

# Check that binary and config both exist
if [ ! -e "$SIMFORAGER_PATH" ]; then
    echo "The executable $SIMFORAGER_PATH does not exist. Exiting."
    exit 1
fi

if [ ! -e "$CONFIG_TEMPLATE_FILE" ]; then
    echo "The configuration file $CONFIG_FILE does not exist. Exiting."
    exit 1
fi

# Convert array task ID
# zero-pad the task ID (01–30) - since task id will be 1-N
SUBSTRATE_NUM=$(printf "%02d" $SLURM_ARRAY_TASK_ID)

# Create an output path based on the job name
OUTPUT_DIR=$HOME/simreef/simreef_outputs/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}/${SUBSTRATE_NUM}/
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

SUBSTRATES_DIR=/carc/scratch/projects/melaniem/melaniem2016345/reef_substrates/SubstrateGenerator/out_spa_459_reps/bmps/
SUBSTRATE_PATH=$SUBSTRATES_DIR/spa_459_rep${SUBSTRATE_NUM}.bmp

CONFIG_FILE=$OUTPUT_DIR/simforager.config
cp $CONFIG_TEMPLATE_FILE $CONFIG_FILE

# Substitute the string <SUBSTRATE_FILE_PATH> in the config file template with the substrate path
sed -i "s|<SUBSTRATE_FILE_PATH>|${SUBSTRATE_PATH}|g" "$CONFIG_FILE"

CMD="srun --mpi=pmi2 $SIMFORAGER_PATH --config $CONFIG_FILE in $OUTPUT_DIR using $SUBSTRATE_PATH"
echo Running $CMD

srun --mpi=pmi2 $SIMFORAGER_PATH --config $CONFIG_FILE 
