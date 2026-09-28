#!/bin/bash

# SLURM script to submit and run Strainify v1.3.1 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=strainify_v1.3.1_simulated_metagenomes
#SBATCH --output=logs/strainify_v1.3.1_simulated_metagenomes.out
#SBATCH --error=logs/strainify_v1.3.1_simulated_metagenomes.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -p long

# Load conda env
source /path/to/conda.sh
conda activate strainify_v1.3.1

# Run strainify v1.3.1
strainify \
  --genome_folder reference_genomes \ # contains either 24 or 22 e.coli reference genomes
  --fastq_folder simulayted_metagenomes \ # simulated metagenomes directory
  --outdir simulated_metagenomes_all_refs_results
