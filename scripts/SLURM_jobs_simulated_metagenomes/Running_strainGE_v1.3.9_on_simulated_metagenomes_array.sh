#!/bin/bash

# SLURM script to run StrainGE v1.3.9 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed
# Should be run when all reference genomes have been kmerized and a pangenome database has been built (see Supplementary Methods)
# Ensure logs directory already exists within the current directory before running the SLURM-scheduled script (mkdir -p logs)

#SBATCH --job-name=strainge_v1.3.9_simulated_metagenomes
#SBATCH --output=/path/to/logs/strainge_v1.3.9_simulated_metagenomes_%a.out
#SBATCH --error=/path/to/logs/strainge_v1.3.9_simulated_metagenomes_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-54%9

# Load the conda environment for StrainGE (adjust the path to your conda installation and environment name)
# See https://github.com/broadinstitute/StrainGE for conda enviornment dependencies
source /path/to/conda.sh
conda activate strainge

# ----------- FULL PATHS (EDIT IF NEEDED) ----------
SAMPLE_LIST="/path/to/simulated_metagenomes.txt" # simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1
PANDB="/path/to/pangenome_db.hdf5"
READDIR="/path/to/simulated_metagenomes_directory"
# --------------------------------------------------

# Get sample prefix from file
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
BASENAME=$(basename "${SAMPLE}")

# Resolve read paths from READDIR (Option A)
R1="${READDIR}/${BASENAME}_R1.fastq.gz"
R2="${READDIR}/${BASENAME}_R2.fastq.gz"

# Make per-sample output directory
OUTDIR="/path/to/simulated_metagenome_results/${BASENAME}"
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Debug: show which files will be used
echo "SAMPLE=${SAMPLE} BASENAME=${BASENAME} R1=${R1} R2=${R2}"
ls -lh "${R1}" "${R2}" || true

# 1) KMERIZE SAMPLE(s)
straingst kmerize \
  -k 31 \
  -o "${BASENAME}.hdf5" \
  "${R1}" "${R2}"

# 2) RUN STRAINGST
straingst run \
  --separate-output \
  -o "${BASENAME}" \
  "${PANDB}" \
  "${BASENAME}.hdf5"

echo "Finished ${BASENAME}"

### End of script