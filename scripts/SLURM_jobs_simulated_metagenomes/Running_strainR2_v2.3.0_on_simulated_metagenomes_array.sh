#!/bin/bash

# SLURM script to submit and run StrainR2 v2.3.0 with BBMap v39.26 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# Ensure that genome preprocessing and database creation (StrainR2DB) with PreProcessR has been completed prior to running this script (see Supplementary Methods)
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=StrainR2_v2.3.0_BBMap_v39.26_simulated_metagenomes_array
#SBATCH --output=logs/StrainR2_v2.3.0_BBMap_v39.26_simulated_metagenomes_array_%a.out
#SBATCH --error=logs/StrainR2_v2.3.0_BBMap_v39.26_simulated_metagenomes_array_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-54%9

# Load the conda environment for StrainR2 (adjust the path to your conda installation and environment name)
# See https://github.com/BisanzLab/StrainR2 for conda enviornment dependencies
# Note our workflow used BBMap v39.26
source /path/to/conda.sh
conda activate strainr2

# Get sample prefix
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" simulated_metagenomes.txt) # simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1

# Make parent output directory
mkdir -p strainR2_BBMap_v39.26_simulated_metagenomes_results

# Run StrainR2
StrainR \
  -1 "${SAMPLE}_R1.fastq.gz" \
  -2 "${SAMPLE}_R2.fastq.gz" \
  -r StrainR2DB \
  -o strainR2_BBMap_v39.26_simulated_metagenomes_results/$(basename ${SAMPLE})

### End of script