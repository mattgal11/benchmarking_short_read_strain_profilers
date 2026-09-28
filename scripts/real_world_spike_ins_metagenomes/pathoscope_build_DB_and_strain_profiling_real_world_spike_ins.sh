#####============================================================#####
# This script is split into two bash scripts to be ran separately:
# 1. PathoScope DB building - the definition of what consitutes the
# "DB-builidng" step are definied here.
# 2. PathoScope strain-level profiling
# This can be used on either the 304 or 294-strain reference database.
#####============================================================#####

#####====================================================#####
# 1. Build the PathoScope database
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate pathoscope

# Limit library thread usage for consistent computational benchmarking
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# CONFIG: edit these variables as appropriate
GENOME_DIR="dRep_0.999_cutoff_genomes" # location of 304 reference genome files (see Methods on how to generate)
DB_DIR="304_ecoli_combined_reference_database" # updated database directory
DB_PREFIX="real_world_spike_ins_combined_ecoli_refs" # prefix for PathoScope to use
TIMING_DIR="computational_benchmarking"
LOG_DIR="logs"

# Create output directories if they don't already exist
mkdir -p "$DB_DIR" "$TIMING_DIR" "$LOG_DIR"

# Time the full DB build for PathoScope:
# 1) concatenate all reference FASTAs into one multi-FASTA
# 2) build the Bowtie2 index from that FASTA
/usr/bin/time -v \
    -o "$TIMING_DIR/${DB_PREFIX}_db_build_time_${THREADS}_threads.txt" \
    sh -c "
        cat '$GENOME_DIR'/*.fna > '$DB_DIR/${DB_PREFIX}.fna' && \
        bowtie2-build '$DB_DIR/${DB_PREFIX}.fna' '$DB_DIR/${DB_PREFIX}'
    " \
    >"$LOG_DIR/${DB_PREFIX}_db_build.out" \
    2>"$LOG_DIR/${DB_PREFIX}_db_build.err"
# end of DB-building script




#####====================================================#####
# 2. Perform strain-level profiling with PathoScope using the
# MAP and ID modules on the ecoli-depleted metagenomes.
# Run this separately after DB-step has completed.
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate pathoscope

# Config: edit these variables as appropriate
FASTQ_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate.
INDEX_DIR="304_ecoli_combined_reference_database" # $DB_DIR from above step
INDEX_PREFIX="real_world_spike_ins_combined_ecoli_refs" # Prefix for PathoScope to use
RESULTS_DIR="real_world_spike_ins_results_pathoscope"
TIMING_DIR="computational_benchmarking"
LOG_DIR="logs"
THREADS=4

# Create results, timing, and log directories if they don't exist
mkdir -p "$RESULTS_DIR" "$TIMING_DIR" "$LOG_DIR"

# Run PathoScope on each sample listed in final_metagenomes.txt
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Missing FASTQ pair for sample: $SAMPLE" >&2
        continue
    fi

    MAP_DIR="$RESULTS_DIR/$SAMPLE/mapOut"
    ID_DIR="$RESULTS_DIR/$SAMPLE/idOut"
    mkdir -p "$MAP_DIR" "$ID_DIR"

    /usr/bin/time -v -o "$TIMING_DIR/${SAMPLE}.time.txt" bash -c "
        pathoscope MAP \
            -numThreads $THREADS \
            -1 '$R1' \
            -2 '$R2' \
            -indexDir '$INDEX_DIR' \
            -targetIndexPrefixes '$INDEX_PREFIX' \
            -outDir '$MAP_DIR' \
            -outAlign '$SAMPLE.sam' \
            -expTag '$SAMPLE' && \
        pathoscope ID \
            -alignFile '$MAP_DIR/$SAMPLE.sam' \
            -fileType sam \
            -outDir '$ID_DIR' \
            -expTag '$SAMPLE'
    " >"$LOG_DIR/${SAMPLE}.out" 2>"$LOG_DIR/${SAMPLE}.err"

done < final_metagenomes.txt # List of samples to process