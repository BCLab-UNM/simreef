#!/bin/bash
#SBATCH --account=2016345
#SBATCH --partition=biocomp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --job-name=spa459
#SBATCH --mail-user=kpowell2@unm.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=92GB
#SBATCH --array=1-30
#SBATCH --output=spa459_%A_%a.out
#SBATCH --error=spa459_%A_%a.err

# -------------------------------
#  Reef Substrate 459 Replicates
# -------------------------------

set -euo pipefail

OUTPUT_DIR=/carc/scratch/projects/melaniem/melaniem2016345/reef_substrates/out_spa_458_reps

# Directories
OUTDIR="${OUTPUT_DIR%/}"
LOGDIR="$OUTDIR/logs"
BMPDIR="$OUTDIR/bmps_5000_karlie"
METADIR="$OUTDIR/meta"

mkdir -p "$LOGDIR" "$BMPDIR" "$METADIR"

# Constants
W=5000
H=5000
PREVALS=("0.25 0.25 0.25 0.25")
SEED_BASE=10000

# Parameters from CSV line 459
ALPHA=1.0
SIGMA=2.0
MIN_PATCH=8000
WARP=22.0874023015592
OCT=4
PERS=0.5199029601650549
LAC=2.207143418394571
FJIT=0.08471110663959631
COMP=0.005695131343509543
BETA=1.0

# Output naming
IDX="${SLURM_ARRAY_TASK_ID}"
SEED=$((SEED_BASE + IDX))

read -r p0 p1 p2 p3 <<< "${PREVALS[0]}"

FNAME="$(printf "%02d.bmp" "$IDX")"
OUTPATH="$BMPDIR/$FNAME"
MANIFEST="$METADIR/manifest_5000_karlie.csv"
LOCKFILE="$METADIR/manifest_5000_karlie.lock"

echo "[$(date)] Starting replicate $IDX with seed $SEED"
echo "[$(date)] OUTPATH: $OUTPATH"

./reef_spa_substrate "$OUTPATH" "$W" "$H" \
  "$ALPHA" "$BETA" "$SIGMA" \
  "$p0" "$p1" "$p2" "$p3" "$SEED" "$MIN_PATCH" \
  "$WARP" "$OCT" "$PERS" "$LAC" "$FJIT" "$COMP"

# Append to manifest with a lock to avoid interleaved lines across array tasks
{
  flock -x 200
  printf '%s\n' "$IDX,$FNAME,$W,$H,$ALPHA,$BETA,$SIGMA,$MIN_PATCH,$WARP,$OCT,$PERS,$LAC,$FJIT,$COMP,$p0,$p1,$p2,$p3,$SEED" >> "$MANIFEST"
} 200>"$LOCKFILE"

echo "[$(date)] Completed replicate $IDX"