#####============================================================#####
# This script is split into two bash scripts to be ran separately:
# 1. StrainGE DB building - the definition of what consitutes the
# "DB-builidng" step are definied here.
# 2. StrainGE strain-level profiling
# This can be used on either the 304 or 294-strain reference database.
#####============================================================#####

#####====================================================#####
# 1. Build the StrainGE database
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate strainge

# Limit library thread usage for consistent benchmarking (StrainGE has no user thread parameter)
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# Create output directories if they don't already exist
mkdir -p computational_benchmarking
mkdir -p dRep_304_ecoli_genomes_hdf5_files

# Build StrainGE database and benchmark the complete workflow (kmerization and database creation)
# The directory "dRep_0.999_cutoff_genomes" should contain the E. coli reference genomes used for benchmarking. 304 reference genomes are provided in the GitHub repo.
/usr/bin/time -v -o computational_benchmarking/strainge_db_build_time_4_threads_export.txt sh -c '
    for f in dRep_0.999_cutoff_genomes/*.fna; do
        base=$(basename "${f%.fna}")

        straingst kmerize \
            -k 31 \
            -o "dRep_304_ecoli_genomes_hdf5_files/${base}.hdf5" \
            "$f"
    done &&

    straingst createdb \
        -o 304_ecoli_genomes_pangenome_db.hdf5 \
        dRep_304_ecoli_genomes_hdf5_files/*.hdf5
'
# end of StrainGE database building



#####====================================================#####
# 2. Strain-level profiling using StrainGE
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate strainge

# Limit library thread usage for consistent benchmarking
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# CONFIG: edit these variables as appropriate
SAMPLE_LIST="final_metagenomes.txt" # List of metagenome sample names.
PANDB="304_ecoli_genomes_pangenome_db.hdf5" # StrainGE database built in the previous step.
READDIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate.
RESULTS_DIR="real_world_spike_ins_strainge_results"
TIMING_DIR="computational_benchmarking"
LOG_DIR="logs"

mkdir -p "$RESULTS_DIR" "$TIMING_DIR" "$LOG_DIR"

while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    R1="${READDIR}/${SAMPLE}_R1.fastq.gz"
    R2="${READDIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Missing FASTQ pair for sample: $SAMPLE" >&2
        continue
    fi

    OUTDIR="${RESULTS_DIR}/${SAMPLE}"
    mkdir -p "$OUTDIR"

    /usr/bin/time -v -o "${TIMING_DIR}/${SAMPLE}.time.txt" sh -c "
        straingst kmerize \
            -k 31 \
            -o ${OUTDIR}/${SAMPLE}.hdf5 \
            ${R1} ${R2} &&

        straingst run \
            --separate-output \
            -o ${OUTDIR}/${SAMPLE} \
            ${PANDB} \
            ${OUTDIR}/${SAMPLE}.hdf5
    " \
    > "${LOG_DIR}/${SAMPLE}.out" \
    2> "${LOG_DIR}/${SAMPLE}.err"

done < "$SAMPLE_LIST"