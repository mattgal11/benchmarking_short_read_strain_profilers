#!/bin/bash

# SLURM script to submit and run Strainify v1.1.0 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# config.yaml must have been updated (see Supplementary Methods) and be present within the current working directory
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=strainify_v1.1.0_simulated_metagenomes_array
#SBATCH --output=logs/strainify_v1.1.0_simulated_metagenomes_array_%a.out
#SBATCH --error=logs/strainify_v1.1.0_simulated_metagenomes_array_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -p short

# Simply update the config.yaml file and submit this script from the main Strainify cloned repository directory

# Load the conda environment for Strainify (adjust the path to your conda installation and environment name)
# See https://github.com/treangenlab/Strainify for conda enviornment dependencies
source /path/to/conda.sh
conda activate strainify

# Run Snakemake and submit config.yaml
snakemake --cores 12 --configfile config.yaml

### End of script