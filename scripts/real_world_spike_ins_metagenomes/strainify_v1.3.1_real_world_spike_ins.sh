#####============================================================#####
# This script runs all ecoli-depleted metagenomes through Strainify
# v1.3.1 for strain-level profiling. The first sample is used to
# generate the precomputed files, which are then reused for all 
# subsequent samples. This is definied as the DB-building step.
# This script can be used on the 304 or 294-strain reference database.
#####============================================================#####

#!/usr/bin/env bash

# Load conda environment
source /path/to/conda.sh
conda activate strainify_v1.3.1

# Limit library thread usage for consistent benchmarking
THREADS=4
export OMP_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS"
export MKL_NUM_THREADS="$THREADS"
export NUMEXPR_NUM_THREADS="$THREADS"
export VECLIB_MAXIMUM_THREADS="$THREADS"

# Config: edit these variables as appropriate
GENOME_DIR="dRep_0.999_cutoff_genomes" # location of 304 reference genome files (see Methods on how to generate, or GitHub repo for identities)
ALL_FASTQ_DIR="combined_SRP373424_ERP005534_final_metagenomes_ecoli_removed" # Directory containing the ecoli-depleted metagenomes. See Methods on how to generate.
OUTROOT="real_world_spike_ins_strainify_v1.3.1_results"
TIMINGDIR="computational_benchmarking"
TMP_FASTQ_DIR="${OUTROOT}/tmp_fastq_dirs"

# Create output directories if they don't already exist
mkdir -p "$OUTROOT" "$TIMINGDIR" "$TMP_FASTQ_DIR"

# Run Strainify on each sample listed in the ALL_FASTQ_DIR
# The first sample will be used to generate the precomputed files.
# This will then be reused for all subsequent samples via use of the --use_precomputed_variants and --precomputed_dir flags.
PRECOMP_DIR=""
FIRST_SAMPLE=1

for r1 in "$ALL_FASTQ_DIR"/*_R1.fastq.gz; do
    sample=$(basename "$r1" | sed -E 's/_R1.*//')
    r2=$(echo "$r1" | sed -E 's/_R1/_R2/')

    sample_dir="$TMP_FASTQ_DIR/$sample"
    mkdir -p "$sample_dir"
    ln -sf "$(realpath "$r1")" "$sample_dir/$(basename "$r1")"
    ln -sf "$(realpath "$r2")" "$sample_dir/$(basename "$r2")"

    if [[ "$FIRST_SAMPLE" -eq 1 ]]; then
        /usr/bin/time -v -o "$TIMINGDIR/${sample}.time.txt" \
            strainify \
                --genome_folder "$GENOME_DIR" \
                --fastq_folder "$sample_dir" \
                --outdir "$OUTROOT/$sample" \
                --max_cpus "$THREADS"

        # Check for the 3 required precomputed files in the output directory to allow progression to the next sample
        M1="$OUTROOT/$sample/filtered_variant_matrix.csv"
        M2="$OUTROOT/$sample/sites.txt"
        M3="$OUTROOT/$sample/reference.fna"

        if [[ ! -f "$M1" || ! -f "$M2" || ! -f "$M3" ]]; then
            echo "ERROR: Strainify did not generate all required precomputed files for $sample in $OUTROOT/$sample." >&2
            echo "Missing status:" >&2
            [[ ! -f "$M1" ]] && echo "  - Missing: filtered_variant_matrix.csv" >&2
            [[ ! -f "$M2" ]] && echo "  - Missing: sites.txt" >&2
            [[ ! -f "$M3" ]] && echo "  - Missing: reference.fna" >&2
            exit 1
        fi

        PRECOMP_DIR="$OUTROOT/$sample"
        FIRST_SAMPLE=0
    else
        /usr/bin/time -v -o "$TIMINGDIR/${sample}.time.txt" \
            strainify \
                --genome_folder "$GENOME_DIR" \
                --fastq_folder "$sample_dir" \
                --use_precomputed_variants \
                --precomputed_dir "$PRECOMP_DIR" \
                --outdir "$OUTROOT/$sample" \
                --max_cpus "$THREADS"
    fi
done

# Clean up temporary symlink directory after all samples finish
rm -rf "$TMP_FASTQ_DIR"

# End of script