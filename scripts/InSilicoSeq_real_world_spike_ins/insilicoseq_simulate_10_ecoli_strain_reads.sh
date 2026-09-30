#!/bin/bash
#SBATCH --job-name=10_strain_ecoli_spikeins_real_world_spike_ins_benchmarking_dataset
#SBATCH --output=logs/10_strain_ecoli_spikeins_real_world_spike_ins_benchmarking_dataset_%A_%a.out
#SBATCH --error=logs/10_strain_ecoli_spikeins_real_world_spike_ins_benchmarking_dataset_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -p short
#SBATCH --array=1-42%11

# -----------------------------------------------------------------------------
# Load conda environment for InSilicoSeq
# -----------------------------------------------------------------------------
source /path/to/conda.sh
conda activate insilicoseq
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# CONFIG
# INPUT_COUNTS should be a TSV with a header and these columns:
# sample_id<TAB>paired_reads
#
# Example:
# sample_id    paired_reads
# sample1      33000000
# sample2      41000000
#
# This corresponds to the total number of paired reads in the metagenome samples of interest (use seqkit stats -a *.gz).
# An example of this file is available on GitHub (SRP373424_ERP005534_total_paired_reads.tsv). 
#
# GENOMES_DIR should contain the 10 E. coli genomes to spike-in, in FASTA format. The accessions numbers (for downloading the genomes from NCBI) are listed below:
# GCA_000025745.1_ASM2574v1_genomic.fna
# GCA_000493595.1_EcoPMV1_genomic.fna
# GCA_001280345.1_ASM128034v1_genomic.fna
# GCA_001900395.1_ASM190039v1_genomic.fna
# GCA_002863845.1_ASM286384v1_genomic.fna
# GCA_002996945.1_ASM299694v1_genomic.fna
# GCA_002812545.1_ASM281254v1_genomic.fna
# GCA_003073815.1_ASM307381v1_genomic.fna
# GCA_000245515.1_ASM24551v1_genomic.fna
# GCA_002948655.1_ASM294865v1_genomic.fna
#
#See Table S8 in the supplementary materials of the paper for more information on these strains.
#
# This script will then generate a 1% spike-in of paired reads from the E. coli genomes specified in GENOMES_DIR, using the abundance profile in ABUND_FILE.
# -----------------------------------------------------------------------------
INPUT_COUNTS="/path/to/SRP373424_ERP005534_total_paired_reads.tsv" # Available on GitHub (contains the total number of paired reads in the metagenome samples of interest)
GENOMES_DIR="/path/to/final_10_strains_ecoli_genomes" # Requires downloading the 10 E. coli genomes from NCBI (see above for accessions)
ABUND_FILE="/path/to/abundance_split_by_contig.txt" # Available on GitHub (InSilicoSeq works on contigs, rather than whole genomes)
OUTDIR="/path/to/SRP373424_ERP005534_ecoli_spike_in_reads"
ISS_MODEL="hiseq"   # Should match the sequencing platform of the metagenome samples
SEED_OFFSET=1000
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Make output and logs directory if they don't exist
# Export thread variables for InSilicoSeq
# -----------------------------------------------------------------------------
mkdir -p logs
mkdir -p "$OUTDIR"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Read sample_id and paired_reads from the TSV
# Assumes one header line, so task 1 corresponds to line 2.
# -----------------------------------------------------------------------------
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))

read -r SAMPLE_ID PAIRED_READS < <(
    awk -F'\t' -v n="$LINE_NUM" 'NR==n {print $1, $2}' "$INPUT_COUNTS"
)

if [[ -z "${SAMPLE_ID:-}" || -z "${PAIRED_READS:-}" ]]; then
    echo "ERROR: Could not read sample_id and paired_reads from line ${LINE_NUM} of ${INPUT_COUNTS}" >&2
    exit 1
fi
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 1% spike-in of the FINAL metagenome
# current single reads = PAIRED_READS * 2 (this represents 99% of final)
# spike-in (1%) = (PAIRED_READS * 2) / 99
# -----------------------------------------------------------------------------
SPIKE_SINGLE=$(( PAIRED_READS * 2 / 99 ))

SEED=$((SEED_OFFSET + SLURM_ARRAY_TASK_ID))
OUT_PREFIX="${OUTDIR}/${SAMPLE_ID}_ecoli_spikein"

echo "Sample:          ${SAMPLE_ID}"
echo "Paired reads:    ${PAIRED_READS}"
echo "Spike singles:   ${SPIKE_SINGLE}"
echo "Seed:            ${SEED}"
echo "Output prefix:   ${OUT_PREFIX}"
echo "Genomes dir:     ${GENOMES_DIR}"
echo "Abundance file:  ${ABUND_FILE}"
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Run InSilicoSeq
# -----------------------------------------------------------------------------
iss generate \
    --model "${ISS_MODEL}" \
    --genomes "${GENOMES_DIR}"/*.fna \
    --abundance_file "${ABUND_FILE}" \
    --n_reads "${SPIKE_SINGLE}" \
    --seed "${SEED}" \
    --cpus "${SLURM_CPUS_PER_TASK:-4}" \
    --debug \
    --compress \
    --output "${OUT_PREFIX}"

echo "Finished ${SAMPLE_ID} at $(date)"
# -----------------------------------------------------------------------------

# End of script