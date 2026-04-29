#!/bin/bash

# This SLURM-scheduled script will simulate the control metagenomes (no K12-MG1655 or O157:H7 str. Sakai) using InSilicoSeq v2.0.1, in triplicate
# See Methods for further information.
# Table S4 contains normalised per-contig relative abundances for species assemblies used to construct the baseline gut microbiome profile.
# Table S5 contains the final InSilicoSeq setup for metagenome simulation, including strain abundances and depth of coverage.
# See Data Summary for .txt --abundance_file (control = baseline community with target strains removed (see Methods))

#SBATCH --job-name=InSilicoSeq_simulated_metagenomes_control
#SBATCH --output=logs/iss_control_%A_%a.out
#SBATCH --error=logs/iss_control_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -p short
#SBATCH --array=1-9%9   # 3 depths * 3 reps = 9 tasks

# Load conda environment for InSilicoSeq (adjust the path to your conda installation and environment name)
# If not set up yet follow conda instructions available on https://github.com/HadrienG/InSilicoSeq
source /path/to/conda.sh
conda activate insilicoseq

# ---------------- CONFIG ------------------------------------------------------------------------------------------------------------------------
REPS=3  # rep1..rep3 (seeds 1,2,3)
N_READS_ARRAY=(40000000 100000000 200000000)  # depths expressed as total single reads (InSilicoSeq --n_reads)
PAIR_LABELS=(20 50 100)   # paired-read million-labels corresponding to N_READS_ARRAY (40M -> 20M pairs -> label 20)
GENOMES_DIR="/path/to/genomes_directory" # should contain assemblies detailed in Table S4, excluding K12-MG1655 and O157:H7 Sakai E. coli strains
OUTDIR="/path/to/simulated_metagenomes" # where we will store our simulated metagenomes that are generated
ABUND_PATH="abundance_files/control.txt"   # single control abundance file (provided; see Data Summary)
# ------------------------------------------------------------------------------------------------------------------------------------------------

# Make directories if they do not already exist
mkdir -p logs
mkdir -p "$OUTDIR"

# Safety: total tasks = depths * reps
ND=${#N_READS_ARRAY[@]}
TOTAL=$(( ND * REPS ))

#Safety check block; ensures there are 9 tasks to be performed, in line with the SLURM array at the start
if [ "$SLURM_ARRAY_TASK_ID" -lt 1 ] || [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL" ]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) out of range (1..$TOTAL)"
  exit 1
fi

# zero-based index
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))

# Map task -> depth index and replicate:
# Group by depth so tasks 1..3 -> depth0 rep1..3, 4..6 -> depth1 rep1..3, etc.
DEP_IDX=$(( IDX / REPS ))          # 0..ND-1
REP=$(( IDX % REPS + 1 ))          # 1..REPS
N_READS=${N_READS_ARRAY[$DEP_IDX]}
PAIR_LABEL=${PAIR_LABELS[$DEP_IDX]}   # 20, 50, or 100 (for filename)
TAG="control"
OUT_PREFIX="${OUTDIR}/HiSeq_${PAIR_LABEL}M_${TAG}_rep${REP}" # output file name structure
SEED=$REP

# Start timestamp (per task)
echo "Started $OUT_PREFIX at $(date '+%Y-%m-%d %H:%M:%S')"

# Print summary
echo "Task $SLURM_ARRAY_TASK_ID -> depth_index=$DEP_IDX N_READS=$N_READS (pairs ${PAIR_LABEL}M) TAG=$TAG REP=$REP SEED=$SEED OUT=${OUT_PREFIX}"
echo "Using genomes dir: $GENOMES_DIR"
echo "Abundance file: $ABUND_PATH"

# Thread count env vars
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

# Run InSilicoSeq using the chosen parameters (see Methods; Table S5)
iss generate \
  --model hiseq \
  --genomes "$GENOMES_DIR"/*.fna \
  --abundance_file "$ABUND_PATH" \
  --n_reads "$N_READS" \
  --seed "$SEED" \
  --cpus "${SLURM_CPUS_PER_TASK:-8}" \
  --debug \
  --compress \
  --output "$OUT_PREFIX"

echo "Finished $OUT_PREFIX at $(date '+%Y-%m-%d %H:%M:%S')"

### End of script