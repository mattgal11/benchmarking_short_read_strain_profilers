# File-driven batch processing loop in Bash to run PanTax v2.0.2 on the simulated metagenomes dataset
# Can be used for full reference database runs or with K12-MG1655 and O157:H7 str. Sakai removed

### Instructions:
# Ensure gurobi_cl for the academic license works before running by setting environment variables (e.g. GUROBI_ROOT, LD_LIBRARY_PATH, and PATH) to the local installation directory
# Ensure you are in the correct directory (where genomes_info.txt and pantax_db are located)
# Ensure the PanTax conda environment has been activated already (conda activate PanTax; see https://github.com/LuoGroup2023/PanTax for dependencies)
# Ensure database construction has been completed (see Supplementary Methods: e.g. pantax -f genomes_info.txt -db --create)

# Create output directories for results and logs
mkdir -p pantax_v.2.0.2_simulated_metagenomes_dataset
mkdir -p pantax_v.2.0.2_simulated_metagenomes_dataset/logs

while IFS= read -r SAMPLE; do
  SAMPLE_ID=$(basename "${SAMPLE}")

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running ${SAMPLE}"

  pantax -f genomes_info.txt -s -p \
    -r "${SAMPLE}_R1.fastq.gz" \
    -r "${SAMPLE}_R2.fastq.gz" \
    -db pantax_db --species --strain \
    --output "pantax_v.2.0.2_simulated_metagenomes_dataset/${SAMPLE_ID}" \
    > "pantax_v.2.0.2_simulated_metagenomes_dataset/logs/${SAMPLE_ID}.log" 2>&1

done < simulated_metagenomes.txt # simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1

# Extract predicted abundance results per run:
echo -e "sampleID\tgenome_ID\tpredicted_abundance" > all_samples_abundance.txt

find pantax_v.2.0.2_simulated_metagenomes_dataset -type f -name '*_strain_abundance.txt' | sort | while IFS= read -r f; do
    sample=$(basename "${f}" _strain_abundance.txt)
    awk -v s="$sample" 'NR>1 {print s"\t"$3"\t"$5}' "$f" >> all_samples_abundance.txt
done

### End of script
