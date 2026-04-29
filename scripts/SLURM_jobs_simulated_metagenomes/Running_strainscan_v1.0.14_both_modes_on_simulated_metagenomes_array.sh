#!/bin/bash

# SLURM script to submit and run StrainScan v1.0.14 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# Ensure that a custom database (reference_genome_database) has been built using strainscan_build prior to running this script (see Supplementary Methods)
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=StrainScan_v1.0.14_simulated_metagenomes_array
#SBATCH --output=logs/StrainScan_v1.0.14_simulated_metagenomes_%a.out
#SBATCH --error=logs/StrainScan_v1.0.14_simulated_metagenomes_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-54%9

# Load the conda environment for StrainScan (adjust the path to your conda installation and environment name)
# See https://github.com/liaoherui/StrainScan for dependencies
source /path/to/conda.sh
conda activate strainscan

# Make parent directory
mkdir -p strainscan_v1.0.14_simulated_metagenomes_results

# Get sample prefix
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" simulated_metagenomes.txt) # simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1

# Run StrainScan to identify bacterial strains in short-reads (paired-end)
strainscan \
  -i "${SAMPLE}_R1.fastq.gz" \
  -j "${SAMPLE}_R2.fastq.gz" \
  -d reference_genome_database \
  -o "strainscan_v1.0.14_simulated_metagenomes_results/$(basename "${SAMPLE}")"

# Alternatively uncomment below and comment above and run StrainScan in low depth mode with probabilistic strain detection on the same samples
# mkdir -p strainscan_v1.0.14_super_low_depth_mode_simulated_metagenomes_results
#strainscan \
#  -i "${SAMPLE}_R1.fastq.gz" \
#  -j "${SAMPLE}_R2.fastq.gz" \
#  -d reference_genome_database \
#  -l 2 \
#  -b 1 \
#  -o "strainscan_v1.0.14_super_low_depth_mode_simulated_metagenomes_results/$(basename "${SAMPLE}")"

### End of script