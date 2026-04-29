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
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running ${SAMPLE}"
  pantax -f genomes_info.txt -s -p \
    -r "${SAMPLE}_R1.fastq.gz" \
    -r "${SAMPLE}_R2.fastq.gz" \
    -db pantax_db --species --strain \
    --output "pantax_v.2.0.2_simulated_metagenomes_dataset/${SAMPLE}" \
    > "pantax_v.2.0.2_simulated_metagenomes_dataset/logs/${SAMPLE}.log" 2>&1
done < simulated_metagenomes.txt
# ^^simulated_metagenomes.txt (provided - see Data Summary) should contain one entry per line, including the relative path to the sample (without _R1/_R2 suffix), e.g., final_metagenomes/HiSeq_20M_K12_50_Sakai_50_rep1

# To extract the predicted abundance results run the following in the results directory:
echo -e "sampleID\tgenome_ID\tpredicted_abundance" > all_samples_abundance.txt

for f in *_strain_abundance.txt; do
    sample=${f%_strain_abundance.txt}
    awk -v s="$sample" 'NR>1 {print s"\t"$3"\t"$5}' "$f" >> all_samples_abundance.txt
done

### End of script