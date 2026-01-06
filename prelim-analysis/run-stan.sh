#!/bin/bash

#SBATCH --output log/st-hawkes-stan.out-%A_%a
#SBATCH --array=1-48
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=48:00:00


# Change directory
cd ..

# Load modules
module purge
module load gdal/3.9.0 proj/9.4.0
module load R

stdbuf -i0 -o0 -e0 command

# Run R script
Rscript prelim-analysis/test_discretetime_st_stan_adm_only_icar_zeroinf.R 
