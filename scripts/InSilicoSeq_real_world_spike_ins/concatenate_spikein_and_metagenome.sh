#!/bin/bash
#SBATCH --job-name=concatenate_spikeins_into_cleaned_metagenomes_SRP373424_ERP005534
#SBATCH --output=logs/concatenate_spikeins_into_cleaned_metagenomes_SRP373424_ERP005534_%A_%a.out
#SBATCH --error=logs/concatenate_spikeins_into_cleaned_metagenomes_SRP373424_ERP005534_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH -p short
#SBATCH --array=1-42%11   # (Adjust to match the number of samples in final_metagenomes.txt)

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
SAMPLE_LIST="/path/to/final_metagenomes.txt" # Available on GitHub (contains the list of metagenome sample names to process)
CLEAN_DIR="/path/to/ERP005534_metagenomes_ecoli_removed" # Directory containing the cleaned metagenome files (see Methods on how to generate)
SPIKE_DIR="/path/to/SRP373424_ERP005534_ecoli_spike_in_reads" # Generated using the insilicoseq_simulate_10_ecoli_strain_reads.sh script available on GitHub
FINAL_DIR="/path/to/SRP373424_ERP005534_final_spiked_metagenomes"

mkdir -p "$FINAL_DIR"
mkdir -p logs
# -----------------------------------------------------------------------------

# Extract the exact sample name for this task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")

echo "Processing Sample: ${SAMPLE}"

# Define explicit paths for this specific sample
CLEAN_R1="${CLEAN_DIR}/${SAMPLE}_ecoli_reads_removed_R1.fastq.gz"
CLEAN_R2="${CLEAN_DIR}/${SAMPLE}_ecoli_reads_removed_R2.fastq.gz"
    
SPIKE_R1="${SPIKE_DIR}/${SAMPLE}_ecoli_spikein_R1.fastq.gz"
SPIKE_R2="${SPIKE_DIR}/${SAMPLE}_ecoli_spikein_R2.fastq.gz"
    
FINAL_R1="${FINAL_DIR}/${SAMPLE}_final_spiked_R1.fastq.gz"
FINAL_R2="${FINAL_DIR}/${SAMPLE}_final_spiked_R2.fastq.gz"

# Safety Check: Ensure all inputs exist before running cat
if [[ -f "$CLEAN_R1" && -f "$SPIKE_R1" ]]; then
    # Concatenate Forward Reads (R1)
    cat "$CLEAN_R1" "$SPIKE_R1" > "$FINAL_R1"
    
    # Concatenate Reverse Reads (R2)
    cat "$CLEAN_R2" "$SPIKE_R2" > "$FINAL_R2"
    
    echo "   SUCCESS: Created spiked files for ${SAMPLE}"
else
    echo "   ERROR: Missing input files for ${SAMPLE}. Skipping!" >&2
fi

echo "Sample ${SAMPLE} processed at $(date)"

# End of script