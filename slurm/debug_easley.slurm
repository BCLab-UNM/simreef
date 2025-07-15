#!/bin/bash
#SBATCH --account 2016345
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --job-name=simreef
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=128GB
#SBATCH --output=debug_%j.out
#SBATCH --error=debug_%j.out

# Setup environment
export UPCXX_SHARED_HEAP_SIZE='16 GB'
export UPCXX_THREADMODE=seq
export UPCXX_CODEMODE=opt
export UPCXX_NETWORK=ibv  # Use 'mpi' if IBV causes issues

# GASNet settings (adjusted for 128 GB RAM)
export GASNET_MAX_SEGSIZE='16 GB'
export GASNET_PHYSMEM_MAX='64 GB'
export GASNET_IBV_PHYSMEM_MAX='64 GB'
export GASNET_USE_ODP=1  # Enable On-Demand Paging for large-memory systems

# Print backtrace on crash
export GASNET_BACKTRACE=1

# Specify the HCA port name on Easley
export GASNET_IBV_PORTS=mlx5_0:1

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp

module load boost/1.88.0-65dn
module load openmpi/4.1.7-e7k3
module load upcxx/2023.9.0-27k4
module load cmake/3.31.6-auvg
module load opencv/4.10.0-ki66
module load libdeflate/1.14-2pby

# Opencv spack package didnt set the LD path
#export LD_LIBRARY_PATH=/opt/spack/opt/spack/linux-sapphirerapids/opencv-4.10.0-ki66vuchmif7c27w2u4hnb3gfynkbbey/lib64:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$OPENCV_LIB64:$LIBDEFLATE_LIB$LD_LIBRARY_PATH

SIMFORAGER_PATH=$HOME/simreef/simforager/install/bin/simforager
CONFIG_FILE=$HOME/simreef/simforager/configs/simforager_reef.config

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
