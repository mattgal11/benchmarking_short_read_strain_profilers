#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------
# -- CONFIG (edit PARENT_DIR and OUTFILE as needed) --
# ----------------------------------------------------
PARENT_DIR="strainGE_results"
OUTFILE="strainGE_relative_strain_abundances.tsv"

# Temporary file for sortable list
TMPFILE=$(mktemp /tmp/strainge_sort.XXXXXX)
trap 'rm -f "$TMPFILE"' EXIT

# Header for final output
echo -e "sample\tstrain\trapct\trel_pct" > "$OUTFILE"

# ----------------------------------------------------
# Build sortable file list (depth → K12 → Sakai)
# ----------------------------------------------------
while IFS= read -r -d '' f; do
  base=$(basename "$f")

  # Extract depth (e.g. 20 from 20M)
  if [[ $base =~ _([0-9]+)M_ ]]; then
    depth_num=${BASH_REMATCH[1]}
  else
    depth_num=9999
  fi

  # Extract K12 proportion
  if [[ $base =~ _K12[_-]?([0-9]+) ]]; then
    k12_num=${BASH_REMATCH[1]}
  else
    k12_num=9999
  fi

  # Extract Sakai proportion
  if [[ $base =~ _Sakai[_-]?([0-9]+) ]]; then
    sakai_num=${BASH_REMATCH[1]}
  else
    sakai_num=9999
  fi

  printf '%s\t%s\t%s\t%s\n' "$depth_num" "$k12_num" "$sakai_num" "$f" >> "$TMPFILE"

done < <(find "$PARENT_DIR" -type f -name '*.strains.tsv' -print0)

# ----------------------------------------------------
# Process files in sorted order
# ----------------------------------------------------
sort -t$'\t' -k1,1n -k2,2n -k3,3n "$TMPFILE" \
| cut -f4- \
| while IFS= read -r file; do

  sample=$(basename "$file" .strains.tsv)

  # ------------------------------------------------
  # If file is empty (e.g. control), write zero row
  # ------------------------------------------------
  if [[ ! -s "$file" ]]; then
    printf "%s\tNA\t0\t0\n" "$sample" >> "$OUTFILE"
    continue
  fi

  # ------------------------------------------------
  # Otherwise compute inter-strain relative abundance
  # ------------------------------------------------
  awk -F'\t' -v sample="$sample" '
  NR==1 {
    for (i=1;i<=NF;i++) {
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
    sum = 0
    for (i=1;i<=n;i++) sum += vals[i]
    for (i=1;i<=n;i++) {
      rel = (sum == 0 ? 0 : (vals[i] / sum) * 100)
      printf "%s\t%s\t%g\t%.2f\n", sample, strains[i], vals[i], rel
    }
  }
  ' "$file" >> "$OUTFILE"

done

echo "Done. Output written to $OUTFILE"