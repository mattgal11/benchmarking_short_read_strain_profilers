#!/bin/bash

# This SLURM-scheduled script will simulate the 5 abundance ratios metagenomes using InSilicoSeq v2.0.1, across 3 depths and in triplicate
# See Methods for further information.
# Table S4 contains normalised per-contig relative abundances for 98 species assemblies used to construct the baseline gut microbiome profile.
# Table S5 contains the final InSilicoSeq setup for metagenome simulation, including strain abundances and depth of coverage.
# See Data Summary for .txt --abundance_file(s)

#SBATCH --job-name=InSilicoSeq_simulated_metagenomes
#SBATCH --output=logs/iss_%A_%a.out
#SBATCH --error=logs/iss_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-45%9   # 5 tags * 3 reps * 3 depths = 45 tasks (safety check later for this) - can adjust %9 if required

# Load conda environment for InSilicoSeq (adjust the path to your conda installation and environment name)
# If not set up yet follow conda instructions available on https://github.com/HadrienG/InSilicoSeq
source /path/to/conda.sh
conda activate insilicoseq

# ---------------- CONFIG ---------------------------------------------------------------------------------------------------------
LIST="abundance_files/sbatch_submission_list.txt"   # this file has one abundance filename per line, 5 lines, e.g. K12_50_Sakai_50
REPS=3  # rep1, rep2,.rep3 (seeds 1,2,3; used for reproducibility)
N_READS_ARRAY=(40000000 100000000 200000000)  # depths expressed as total single reads (InSilicoSeq --n_reads)
PAIR_LABELS=(20 50 100)   # paired-read million-labels corresponding to N_READS_ARRAY (40M -> 20M pairs -> label 20)
GENOMES_DIR="/path/to/genomes_directory"  # genomes_directory should contain the assemblies detailed in Table S4
OUTDIR="/path/to/simulated_metagenomes" # where we will store our simulated metagenomes that are generated
# ---------------------------------------------------------------------------------------------------------------------------------

# Make directories if they do not already exist
mkdir -p logs
mkdir -p "$OUTDIR"

# Computes TOTAL=abundances × replicates × depths (or 5 x 3 x 3, for 45 files total files)
NLINES=$(wc -l < "$LIST")
ND=${#N_READS_ARRAY[@]}
TOTAL=$(( NLINES * REPS * ND ))

#Safety check block; ensures there are 45 tasks to be performed, in line with the SLURM array at the start
if [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) out of range (1..$TOTAL)"
  exit 1
fi

# zero-based index (SLURM array IDs start at 1 but programming logic is easier with 0-based indexing for bash)
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))

# determine depth index, then remainder for line & replicate
BLOCK=$(( NLINES * REPS ))          # number of tasks per depth (e.g., 1..15 are depth 0 (20M))
DEP_IDX=$(( IDX / BLOCK ))          # 0..ND-1
REM=$(( IDX % BLOCK ))              # remainder within the chosen depth block (0..BLOCK-1).
LINE=$(( REM % NLINES + 1 ))        # 1..NLINES (from sbatch_submission_list.txt)
REP=$(( REM / NLINES + 1 ))         # 1..REPS
TAG=$(sed -n "${LINE}p" "$LIST")               # prints the LINE-th line of list.txt
ABUND_PATH="abundance_files/${TAG}.txt"        # constructs the abundance filename by appending .txt in the abundance_files/ folder.
N_READS=${N_READS_ARRAY[$DEP_IDX]}
PAIR_LABEL=${PAIR_LABELS[$DEP_IDX]}   # 20, 50, or 100 (for filename)
OUT_PREFIX="${OUTDIR}/HiSeq_${PAIR_LABEL}M_${TAG}_rep${REP}"       # output prefix iss command uses; embeds depth label, tag, and replicate - e.g., HiSeq_20M_K12_50_Sakai_50_rep1
SEED=$REP                              # sets RNG seed to 1,2,3

# Print summary log to verify what each task is doing
echo "Task $SLURM_ARRAY_TASK_ID -> depth_index=$DEP_IDX N_READS=$N_READS (pairs ${PAIR_LABEL}M) TAG=$TAG REP=$REP SEED=$SEED OUT=${OUT_PREFIX}"
echo "Using genomes dir: $GENOMES_DIR"
echo "Abundance file: $ABUND_PATH"

# Thread count enviornment variables are set so that ISS libraries respect the CPU allocation provided by SLURM
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Run InSilicoSeq using the chosen parameters (see Methods; Table S5)
iss generate \
  --model hiseq \
  --genomes $GENOMES_DIR/*.fna \
  --abundance_file "$ABUND_PATH" \
  --n_reads "$N_READS" \
  --seed "$SEED" \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --debug \
  --compress \
  --output "$OUT_PREFIX"

echo "Finished $OUT_PREFIX at $(date)"

### End of script