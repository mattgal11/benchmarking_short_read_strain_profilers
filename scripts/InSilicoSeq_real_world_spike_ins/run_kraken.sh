#!/bin/bash
#SBATCH --job-name=kraken2_ecoli_free_metagenomes_SRP373424_ERP005534
#SBATCH --output=/logs/kraken2_ecoli_free_metagenomes_SRP373424_ERP005534_%A_%a.out
#SBATCH --error=/logs/kraken2_ecoli_free_metagenomes_SRP373424_ERP005534_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-42%11   # <-- Adjust to match the number of samples in final_metagenomes.txt

# -----------------------------------------------------------------------------
# Load conda environment for Kraken2
# -----------------------------------------------------------------------------
source /path/to/conda.sh
conda activate kraken2
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CONFIG: Adjust the file paths here for your given kraken2 run.
# -----------------------------------------------------------------------------
KRAKEN2_DB="/path/to/kraken2_db" # Available from https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz
METAGENOME_DIR="/path/to/SRP373424_and_ERP005534_metagenomes" # Directory contains either the original or the cleaned metagenome files
SAMPLE_LIST="/path/to/final_metagenomes.txt" # Available on GitHub (contains the list of metagenome sample names to process)
OUTDIR="/path/to/kraken2_SRP373424_ERP005534_results"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Make output and logs directory if they don't exist
# -----------------------------------------------------------------------------
mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/logs"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Get sample prefix from sample list file based on SLURM_ARRAY_TASK_ID
# Resolve sample read paths from METAGENOME_DIR
# Show which samples are being processed for debugging purposes
# -----------------------------------------------------------------------------
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

R1="${METAGENOME_DIR}/${SAMPLE}_1.fastq.gz"
R2="${METAGENOME_DIR}/${SAMPLE}_2.fastq.gz"

echo "SAMPLE=${SAMPLE} R1=${R1} R2=${R2} with ${KRAKEN2_DB} database"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run Kraken2 on each of the metagenome samples listed in SAMPLE_LIST
# -----------------------------------------------------------------------------
kraken2 \
    --db  "${KRAKEN2_DB}" \
    --paired "${R1}" "${R2}" \
    --output "${OUTDIR}/${SAMPLE}.kraken" \
    --report "${OUTDIR}/${SAMPLE}.kreport" \
    --use-names \
    --threads 12

echo "Finished ${SAMPLE} at $(date)"
# -----------------------------------------------------------------------------

# End of script