#####============================================================#####
# This script is used to extract and normalise the StrainGE results 
# for the real-world spike-in benchmarking dataset.
# The output TSV file contains the sample name, strain ID, absolute 
# abundance (rapct), and relative abundance (rel_pct) for each strain.
#####============================================================#####

#!/usr/bin/env bash
set -euo pipefail

##### CONFIG #####
PARENT_DIR="real_world_spike_ins_strainge_results" # directory containing the StrainGE results for the real-world spike-in benchmarking dataset
OUTFILE="${PARENT_DIR}/strainge_real_world_spike_ins_all_strains_abundances.tsv"
##################

# Header for final output
printf "sample\tstrain\trapct\trel_pct\n" > "$OUTFILE"

# Process all StrainGE result files in alphabetical order
find "$PARENT_DIR" -type f -name '*.strains.tsv' | sort | while IFS= read -r file; do

    sample=$(basename "$file" .strains.tsv)

    # If file is empty, write zero row
    if [[ ! -s "$file" ]]; then
        printf "%s\tNA\t0\t0\n" "$sample" >> "$OUTFILE"
        continue
    fi

    # Otherwise compute inter-strain relative abundance
    awk -F'\t' -v sample="$sample" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            h = $i
            gsub(/^[ \t]+|[ \t]+$/, "", h)
            if (h == "strain") s_col = i
            if (h == "rapct")  r_col = i
        }
        next
    }
    {
        strain = $s_col
        rapct = ($r_col + 0)
        strains[++n] = strain
        vals[n] = rapct
    }
    END {
        if (n == 0) {
            printf "%s\tNA\t0\t0\n", sample
        } else {
            sum = 0
            for (i=1; i<=n; i++) sum += vals[i]
            for (i=1; i<=n; i++) {
                rel = (sum == 0 ? 0 : (vals[i] / sum) * 100)
                printf "%s\t%s\t%g\t%.2f\n", sample, strains[i], vals[i], rel
            }
        }
    }
    ' "$file" >> "$OUTFILE"

done

echo "Done. Output written to $OUTFILE"