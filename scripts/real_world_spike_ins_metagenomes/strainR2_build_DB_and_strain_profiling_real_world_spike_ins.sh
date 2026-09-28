#####============================================================#####
# This script is split into two bash scripts to be ran separately:
# 1. StrainR2 DB building - the definition of what consitutes the
# "DB-builidng" step are definied here.
# 2. StrainR2 strain-level profiling
# This can be used on either the 304 or 294-strain reference database.
#####============================================================#####

#####====================================================#####
# 1. Build the StrainR2 database
#####====================================================#####

#!/bin/bash

# Load your conda environment
source /path/to/conda.sh
conda activate strainr2

# Limit library thread usage for consistent benchmarking
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# Run PreProcessR (build DB step) with computational benchmarking
# The directory "dRep_0.999_cutoff_genomes" should contain the E. coli reference genomes used for benchmarking. 304 reference genomes are provided in the GitHub repo.
/usr/bin/time -v -o computational_benchmarking/preprocessr_db_build_time.txt \
    PreProcessR -i dRep_0.999_cutoff_genomes \
    -o strainR2_real_world_spike_ins_db
#end of DB building step




#####====================================================#####
# 2. Strain-level profiling with StrainR2
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate strainr2

# Export thread variables
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# CONFIG: edit these variables as appropriate
FASTQ_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate.
DB="strainR2_real_world_spike_ins_db" # DB directory generated in the previous PreProcessR step.
RESULTS_DIR="real_world_spike_ins_strainr2_results"
TIMING_DIR="computational_benchmarking"
LOG_DIR="logs"

# Create output directories if they don't exist
mkdir -p "$RESULTS_DIR" "$TIMING_DIR" "$LOG_DIR"

# Run StrainR2 on each sample listed in final_metagenomes.txt
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Missing FASTQ pair for sample: $SAMPLE" >&2
        continue
    fi

    /usr/bin/time -v -o "$TIMING_DIR/${SAMPLE}.time.txt" \
        StrainR \
            -1 "$R1" \
            -2 "$R2" \
            -r "$DB" \
            -o "$RESULTS_DIR/$SAMPLE" \
            -t "$THREADS" \
            -m 64 \
            > "$LOG_DIR/${SAMPLE}.out" \
            2> "$LOG_DIR/${SAMPLE}.err"
done < final_metagenomes.txt # List of final metagenome samples to process (provided in GitHub repo).
# end of StrainR2 strain-level profiling step

# StrainR2 requires the user to specify the maximum memory available to the underlying BBMap alignment step (-m). 
# The default allocation of 8 GB was insufficient for the benchmark dataset, so this parameter was increased across 12-64 GB to allow successful completion (see Fig.S7).
# Peak memory usage was measured independently using /usr/bin/time -v and is reported as the observed maximum resident set size.