#!/bin/bash
jid=$(sbatch --parsable substrate_generator_array.slurm)
echo "Submitted array job $jid"
dep=$(sbatch --parsable --dependency=afterok:$jid substrate_generator_mosaic.slurm)
echo "Submitted mosaic job $dep (depends on $jid)"
