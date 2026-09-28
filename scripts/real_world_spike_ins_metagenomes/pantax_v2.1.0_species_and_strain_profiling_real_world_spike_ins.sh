#####============================================================#####
# This script is split into two parts:
# 1. PanTax species level profiling with --next flag (SLURM array job)
# 2. PanTax strain level profiling with --next flag (bash scirpt)
#####============================================================#####

#####====================================================#####
# Part 1: PanTax read-only database benchmarking array SLURM
# array job script. Uses --species and --next flags.
# This is run separately for two reasons:
# 1. A higher number of threads is needed for Giraffe mapping 
# step, which is the limiting factor in the overall runtime.
# 2. --strain flag for strain-level profiiling requires gurobi 
# license, which is not available on the compute nodes (no 
# internet access and Gurobi academic license ID is restricted
# to one device - see https://www.gurobi.com/academics).
#####====================================================#####

#!/bin/bash

#SBATCH --job-name=pantax_read_only_DB_array_test
#SBATCH --output=logs/pantax_read_only_DB_array_test_%A_%a.out
#SBATCH --error=logs/pantax_read_only_DB_array_test_%A_%a.err
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH -p long
#SBATCH --array=1-42%12

# ============================================================
# Load conda environment
# ============================================================
source /path/to/conda.sh
conda activate pantax_v2.1.0

# ============================================================
# Configuration
# ============================================================
SAMPLE_LIST="final_metagenomes.txt" # provided in GitHub repo
METAGENOME_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate
DB_DIR="pantax_db" # can be the DB for 304 or 294-strain reference database. See supplementary methods for how to generate the DB.
RESULTS_DIR="pantax_read_only_DB_array_test_results" # output directory for species-level results
TIMING_DIR="computational_benchmarking_read_only_DB_array_test" # output directory for computational benchmarking results
LOG_DIR="${RESULTS_DIR}/logs"
GENOMES_INFO="genomes_info.txt"
THREADS=32

# ============================================================
# Get sample corresponding to this array task
# ============================================================
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

if [[ -z "$SAMPLE" ]]; then
    echo "ERROR: No sample found for array index ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample:      ${SAMPLE}"

# ============================================================
# Paths
# ============================================================
R1="${METAGENOME_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${METAGENOME_DIR}/${SAMPLE}_R2.fastq.gz"

TMP_DIR="${RESULTS_DIR}/${SAMPLE}_tmp"

# ============================================================
# Create output directories
# ============================================================
mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$TIMING_DIR" "$TMP_DIR"

# ============================================================
# Check input FASTQs
# ============================================================
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing FASTQ pair for $SAMPLE" >&2
    echo "R1: $R1" >&2
    echo "R2: $R2" >&2
    exit 1
fi

# ============================================================
# Check database
# We ensure the DB is read-only to avoid accidental overwriting 
# of the large database during benchmarking.
# ============================================================
for index in \
    "$DB_DIR/reference_pangenome.giraffe.gbz" \
    "$DB_DIR/reference_pangenome.min" \
    "$DB_DIR/reference_pangenome.dist"
do
    if [[ ! -r "$index" ]]; then
        echo "ERROR: Index is not readable: $index" >&2
        exit 1
    fi
    if [[ -w "$index" ]]; then
        echo "ERROR: Index is writable: $index" >&2
        exit 1
    fi
done

# ============================================================
# Run PanTax (species-level profiling with --next flag)
# ============================================================
/usr/bin/time -v \
    -o "$TIMING_DIR/${SAMPLE}.time.txt" \
    pantax \
        -f "$GENOMES_INFO" \
        -s \
        -p \
        -r "$R1" \
        -r "$R2" \
        --db "$DB_DIR" \
        --species \
        --next \
        -t "$THREADS" \
        -T "$TMP_DIR" \
        -o "$RESULTS_DIR/$SAMPLE" \
    > "$LOG_DIR/${SAMPLE}.log" \
    2>&1

EXIT_CODE=$?

echo "Sample: $SAMPLE"
echo "PanTax exit code: $EXIT_CODE"

if [[ $EXIT_CODE -eq 0 ]]; then
    if [[ ! -s "$TMP_DIR/gfa_mapped.gaf" ]]; then
        echo "ERROR: PanTax completed but retained GAF is missing/empty" >&2
        exit 2
    fi

    echo "Stage 1 completed successfully"
    echo "Retained tmp directory: $TMP_DIR"
    du -sh "$TMP_DIR"
fi

exit "$EXIT_CODE"



#####====================================================#####
# Part 2: PanTax strain profiling benchmarking bash script
# This script is run on a machine with an internet connection, 
# as the strain profiling step requires pinging Gurobi for the
# license ID, which is also restricted to one device.
# This step is far less computationally intensive and can
# easily be ran using e.g., 4 threads.
#####====================================================#####

#!/bin/bash

# ============================================================
# Load PanTax + export Gurobi environment
# ============================================================

source /path/to/conda.sh
conda activate pantax_v2.1.0

export GUROBI_ROOT=/well/bag/users/wfr051/databases/gurobi1103/linux64
export LD_LIBRARY_PATH="${GUROBI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
export PATH="${GUROBI_ROOT}/bin:${PATH}"

# ============================================================
# Configuration
# ============================================================
SAMPLE_LIST="final_metagenomes.txt"
METAGENOME_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed"
DB_DIR="pantax_db"
STAGE1_DIR="pantax_read_only_DB_array_test_results" # IMPORTANT: This must be the Stage-1 results directory containing SAMPLE_tmp/
RESULTS_DIR="pantax_read_only_DB_strain_resume_final_benchmark_results" # output directory for (final) strain-level results
LOG_DIR="${RESULTS_DIR}/logs"
TIMING_DIR="computational_benchmarking_pantax_strain_resume"
GENOMES_INFO="genomes_info.txt"
THREADS=4

mkdir -p "$RESULTS_DIR" "$LOG_DIR" "$TIMING_DIR"

# ============================================================
# Resume each sample at strain profiling
# ============================================================
while IFS= read -r SAMPLE; do

    [[ -z "$SAMPLE" ]] && continue

    R1="${METAGENOME_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${METAGENOME_DIR}/${SAMPLE}_R2.fastq.gz"

    TMP_DIR="${STAGE1_DIR}/${SAMPLE}_tmp"

    echo "========================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')]"
    echo "Resuming strain profiling: $SAMPLE"
    echo "Tmp: $TMP_DIR"
    echo "========================================"

    # Check Stage-1 state exists
    if [[ ! -d "$TMP_DIR" ]]; then
        echo "ERROR: Missing retained tmp directory for $SAMPLE" >&2
        continue
    fi

    if [[ ! -s "$TMP_DIR/gfa_mapped.gaf" ]]; then
        echo "ERROR: Missing/empty retained GAF for $SAMPLE" >&2
        continue
    fi

    # Check FASTQs still exist
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "ERROR: Missing FASTQs for $SAMPLE" >&2
        continue
    fi

    # ========================================================
    # Stage 2: strain profiling only
    # ========================================================

    /usr/bin/time -v \
        -o "$TIMING_DIR/${SAMPLE}.time.txt" \
        pantax \
            -f "$GENOMES_INFO" \
            -s \
            -p \
            -r "$R1" \
            -r "$R2" \
            --db "$DB_DIR" \
            --strain \
            --next \
            -t "$THREADS" \
            -T "$TMP_DIR" \
            -o "$RESULTS_DIR/$SAMPLE" \
        > "$LOG_DIR/${SAMPLE}.log" \
        2>&1

    EXIT_CODE=$?

    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished $SAMPLE"
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] FAILED $SAMPLE (exit $EXIT_CODE)" >&2
    fi

done < "$SAMPLE_LIST"

echo
echo "All available Stage-1 samples processed."
echo "PanTax finished running on all samples at $(date '+%Y-%m-%d %H:%M:%S')"