#####============================================================#####
# This script is split into two bash scripts to be ran separately:
# 1. StrainScan DB building - the definition of what consitutes the
# "DB-builidng" step are definied here. This can be done using the 
# memory efficient mode (-e 1) (default: 0).
# 2. StrainScan strain-level profiling. This can be using super-low-
# depth mode, as used across our experiments. (-l 2) (default: 0).
# This can be used on either the 304 or 294-strain reference database.
#####============================================================#####

#####====================================================#####
# 1. Build DB with StrainScan
#####====================================================#####

#!/bin/bash

# Load conda
source /path/to/conda.sh
conda activate strainscan

# Limit library thread usage for consistent benchmarking
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# Run StrainScan build DB step with computational benchmarking
# The directory "dRep_0.999_cutoff_genomes" should contain the E. coli reference genomes used for benchmarking. 304 reference genomes are provided in the GitHub repo.
/usr/bin/time -v -o computational_benchmarking/strainscan_db_build_time_dingo_4_threads.txt \
    strainscan_build \
        -i dRep_0.999_cutoff_genomes/ \
        -o real_world_spike_ins_304_ecoli_strain_db/ \
        -t "$THREADS"
        # -e 1 # optional: use memory efficient mode (default: 0)
# end of StrainScan DB building step



#####====================================================#####
# 2. Strain-level profiling with StrainScan
#####====================================================#####

#!/bin/bash

# Load conda environment
source /path/to/conda.sh
conda activate strainscan

# CONFIG: edit these variables as appropriate
FASTQ_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate.
DB="real_world_spike_ins_304_ecoli_strain_db" # StrainScan DB directory built in the previous step
RESULTS_DIR="real_world_spike_ins_results_strainscan"
TIMING_DIR="computational_benchmarking"
LOG_DIR="logs"
THREADS=4

# Export thread variables for StrainScan as it does not have a --threads flag available
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# Make results, timing and log directories if they don't exist
mkdir -p "$RESULTS_DIR" "$TIMING_DIR" "$LOG_DIR"

# Loop through each sample in the input file and run StrainScan
while IFS= read -r SAMPLE; do
    [[ -z "$SAMPLE" ]] && continue

    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Missing FASTQ pair for sample: $SAMPLE" >&2
        continue
    fi

    /usr/bin/time -v -o "$TIMING_DIR/${SAMPLE}.time.txt" \
        strainscan \
            -i "$R1" \
            -j "$R2" \
            -d "$DB" \
            # -l 2 # optional: use super-low-depth mode (default: 0)
            # -b 1 # optional: enable probabalistic strain detection (default: 0)
            -o "$RESULTS_DIR/$SAMPLE" \
        >"$LOG_DIR/${SAMPLE}.out" \
        2>"$LOG_DIR/${SAMPLE}.err"

done < final_metagenomes.txt # List of final metagenome samples to process (provided in GitHub repo).
# End of StrainScan strain-level profiling step.