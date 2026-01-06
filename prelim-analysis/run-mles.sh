#!/bin/bash

#SBATCH --output log/st-hawkes-mle.out-%A_%a
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=48:00:00

# Change directory
cd ..

# Load modules
module purge
module load gdal/3.9.0 proj/9.4.0
module load R

stdbuf -i0 -o0 -e0 command


# Define the array of values
VALUES=("12" "26" "52")

# Pick the value corresponding to this job's array index
VAL=${VALUES[$SLURM_ARRAY_TASK_ID-1]}

echo "Running job $SLURM_ARRAY_TASK_ID with value = $VAL"
Rscript prelim-analysis/test_discretetime_st_mle_adm_only_zeroinf_powerlaw.R "$VAL"
