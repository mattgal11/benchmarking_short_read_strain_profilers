#!/bin/bash

# SLURM script to run PathoScope v2.0.7 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=pathoscope_v2.0.7_simulated_metagenomes_pathoMAP_array
#SBATCH --output=logs/pathoscope_v2.0.7_simulated_metagenomes_pathoMAP_array_%a.out
#SBATCH --error=logs/pathoscope_v2.0.7_simulated_metagenomes_pathoMAP_array_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH -p short
#SBATCH --array=1-54%9 # adjust if you have more or fewer samples (54 in the simulated metagenomes dataset)

# Load the conda environment for PathoScope (adjust the path to your conda installation and environment name)
# See https://github.com/PathoScope/PathoScope for conda enviornment dependencies
source /path/to/conda.sh
conda activate pathoscope

# Read the sample prefix for this array task. simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1
SAMPLE_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" simulated_metagenomes.txt)
# strip whitespace
SAMPLE=$(echo "${SAMPLE_LINE}" | tr -d '\r\n\t ')
# get basename in case SAMPLE contains a path
BASENAME=$(basename "${SAMPLE}")

# paths to input files in the current directory (adjust if files are in a different directory)
R1="${SAMPLE}_R1.fastq.gz"
R2="${SAMPLE}_R2.fastq.gz"

# output directory for this sample
OUTDIR="/path/to/pathoscope_directory/simulated_metagenomes/MapOut/${BASENAME}"
mkdir -p "${OUTDIR}"

# reference (adjust path; assumes the reference files are in the current directory for which a Bowtie2 index has been built; see Supplementary Methods)
TARGET_REFS="simulated_metagenomes_refs.fna"

# run PathoScope MAP
pathoscope MAP \
  -1 "${R1}" \
  -2 "${R2}" \
  -targetRefFiles "${TARGET_REFS}" \
  -outDir "${OUTDIR}" \
  -outAlign "${OUTDIR}/${BASENAME}.sam" \
  -expTag "${BASENAME}"

### End of script