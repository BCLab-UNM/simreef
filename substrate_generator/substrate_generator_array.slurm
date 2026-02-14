#!/bin/bash
#SBATCH --account=2016345
#SBATCH --partition=biocomp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --job-name=substrate_gen
#SBATCH --mail-user=mfricke@unm.edu
#SBATCH --mail-type=FAIL
#SBATCH --mem=2GB
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --array=0-511

module load gcc/11.3.0
cd $SLURM_SUBMIT_DIR

PARAMS_FILE=out_spa_512/meta/params_512.csv

if [[ ! -f $PARAMS_FILE ]]; then
    echo "ERROR: $PARAMS_FILE not found. Run generate_params.py first."
    exit 1
fi

OUT_DIR=out_spa_512/bmps
mkdir -p $OUT_DIR logs

# Pick the correct parameter row
PARAMS=$(awk -F, -v n=$((SLURM_ARRAY_TASK_ID+2)) 'NR==n{print}' $PARAMS_FILE)
IFS=',' read alpha sigma mp warp oct pers lac fjit comp <<< "$PARAMS"

# Constants
p0=0.25; p1=0.25; p2=0.25; p3=0.25
seed=$((12345 + SLURM_ARRAY_TASK_ID))
BETA=1.0
W=5000
H=5000

fname=$(printf "spa_%03d_a%.2f_s%.1f_m%05d_w%02.0f.bmp" $SLURM_ARRAY_TASK_ID $alpha $sigma $mp $warp)
out="$OUT_DIR/$fname"

echo "[$SLURM_ARRAY_TASK_ID] Generating $out"
./reef_spa_substrate "$out" $W $H "$alpha" $BETA "$sigma" \
  $p0 $p1 $p2 $p3 $seed $mp "$warp" "$oct" "$pers" "$lac" "$fjit" "$comp"
