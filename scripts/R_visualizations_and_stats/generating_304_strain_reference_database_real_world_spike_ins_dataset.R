# ============================================================
# Additional R script: generating the 304-strain reference
# database used in the real-world spike-ins benchmarking dataset
#
# This is based on the E. coli reference genomes available in
# the Abram et al., (2021) publication.
# DOI: 10.1038/s42003-020-01626-5
#
# This script includes steps up to point of dereplication, 
# where we are left with 526 complete genomes to be used in
# dRep, using the parameters outlined in Supp. Methods S2.
# ============================================================

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readxl)

# Read in the NCBI datasets assembly stats dataframe for the ecoli reference genomes from Abram et al., (2021) paper
# See Supplementary Methods S2 for further details, or the file below itself (has additional sheet called "generation_info.txt")
abram_assembly_stats <- read_xlsx(
  "abram_et_al_assembly_stats.xlsx", sheet = "abram_et_al_assembly_stats") # <-- Available on FigShare

# Read in the original Abram et al., (2021) ecoli-only publication data for phylogroup information
abram_phylo <- read_xlsx(
  "abram_et_al_2021_ecoli_reference_genomes.xlsx") # <-- This requires downloading from the original publication by Abram et al., (2021). https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-020-01626-5/MediaObjects/42003_2020_1626_MOESM3_ESM.csv

# Merge the two dataframes on the assembly accession column to get a combined dataframe with assembly stats and phylogroup information
abram_assembly_stats_phylo <- abram_assembly_stats %>%
  left_join(abram_phylo %>%
      select(`Assembly Accession`, Phylogroup, `#Organism/Name`, Strain, `Size (Mb)`, Bioproject, Biosample),
      by = c("accession" = "Assembly Accession"))

# Filter the merged table for complete genomes only and remove rows where checkM scores are NA
abram_complete_genomes <- abram_assembly_stats_phylo %>%
  filter(assembly_level == "Complete Genome",
    !is.na(checkm_completeness),
    !is.na(checkm_contamination))


# Plot scatter plot of checkM contamination vs checkM completeness
p_complete_genomes_checkM <- ggplot(abram_complete_genomes, aes(x = checkm_contamination , y = checkm_completeness)) +
  geom_point() +
  scale_x_log10(labels = scales::comma) +
  labs(title = "CheckM contamination (%) vs CheckM Completeness (%) for complete genomes",
       x = "CheckM contamination (%, log10 scale)",
       y = "CheckM Completeness (%)") +
  theme_minimal()

print(p_complete_genomes_checkM)


# Generate retention table across checkM completeness and contamination cutoffs for the complete genomes
checkm_comp <- abram_complete_genomes$checkm_completeness
checkm_cont <- abram_complete_genomes$checkm_contamination

keep <- is.finite(checkm_comp) & is.finite(checkm_cont)
checkm_comp <- checkm_comp[keep]
checkm_cont <- checkm_cont[keep]

comp_cutoffs <- c(95, 97, 98, 99, 99.5)
cont_cutoffs <- c(0.5, 1, 2, 3, 5)

retention_matrix <- outer(comp_cutoffs, cont_cutoffs, Vectorize(function(cmp, cont) {
  sum(checkm_comp >= cmp & checkm_cont <= cont)
}))

rownames(retention_matrix) <- paste0("comp>=", comp_cutoffs)
colnames(retention_matrix) <- paste0("cont<=", cont_cutoffs)

retention_matrix

n_total <- length(checkm_comp) # get the percentage too
retention_percent_matrix <- round(100 * retention_matrix / n_total, 1)
retention_percent_matrix

write.csv(
  retention_percent_matrix, "abram_et_al_complete_genomes_checkM_cut_offs.csv", row.names = FALSE)

  
# Long table version of the above
retention_table <- expand.grid(
  completeness_cutoff = comp_cutoffs,
  contamination_cutoff = cont_cutoffs)

retention_table$retained <- mapply(function(cmp, cont) {
  sum(checkm_comp >= cmp & checkm_cont <= cont)
}, retention_table$completeness_cutoff, retention_table$contamination_cutoff)

retention_table$retained_percent <- round(100 * retention_table$retained / length(checkm_comp), 1)

retention_table
write.csv(
  retention_table, "abram_et_al_complete_genomes_checkM_cut_offs_long_table.csv", row.names = FALSE)


# Heatmap plot of the above retention percent matrix
library(reshape2)

heatmap_df <- melt(retention_percent_matrix)
colnames(heatmap_df) <- c(
  "completeness_cutoff",
  "contamination_cutoff",
  "retained_percent")

p_complete_genomes_retained_heatmap <- ggplot(heatmap_df,
  aes(x = contamination_cutoff, y = completeness_cutoff,fill = retained_percent)) +
  geom_tile() +
  geom_text(aes(label = retained_percent), size = 4) +
  scale_fill_viridis_c(name = "% retained") +
  labs(
    x = "CheckM contamination (%) cutoff",
    y = "CheckM completeness (%) cutoff",
    title = "Genome retention under CheckM QC thresholds") +
  theme_minimal()

print(p_complete_genomes_retained_heatmap)


# do the same at the phylogroup level to check retention per phylogroup
library(dplyr)
library(tidyr)

comp_cutoffs <- c(95, 97, 98, 99, 99.5)
cont_cutoffs <- c(0.5, 1, 2, 3, 5)

  # total genomes per phylogroup
phylo_totals <- abram_complete_genomes %>%
  count(Phylogroup, name = "Total")

  # counts retained for each completeness/contamination cutoff combination
phylo_counts <- expand_grid(
  completeness_cutoff = comp_cutoffs,
  contamination_cutoff = cont_cutoffs
) %>%
  rowwise() %>%
  mutate(
    data = list(
      abram_complete_genomes %>%
        filter(
          checkm_completeness >= completeness_cutoff,
          checkm_contamination <= contamination_cutoff
        ) %>%
        count(Phylogroup, name = "retained")
    )
  ) %>%
  unnest(data) %>%
  ungroup() %>%
  left_join(phylo_totals, by = "Phylogroup") %>%
  mutate(retained_percent = round(100 * retained / Total, 1))

# fill in missing combinations with 0s

phylo_counts <- phylo_counts %>%
  complete(
    Phylogroup,
    completeness_cutoff,
    contamination_cutoff,
    fill = list(
      retained = 0,
      retained_percent = 0
    )
  )

# Wide table of retained counts
phylo_count_table <- phylo_counts %>%
  mutate(cutoff = paste0("comp_", completeness_cutoff, "_cont_", contamination_cutoff)) %>%
  select(Phylogroup, Total, cutoff, retained) %>%
  pivot_wider(
    names_from = cutoff,
    values_from = retained,
    values_fill = 0)

# Wide table of retained percentages
phylo_percent_table <- phylo_counts %>%
  mutate(cutoff = paste0("comp_", completeness_cutoff, "_cont_", contamination_cutoff)) %>%
  select(Phylogroup, Total, cutoff, retained_percent) %>%
  pivot_wider(
    names_from = cutoff,
    values_from = retained_percent,
    values_fill = 0)


# Plot the above as a phylogroup-specific heatmap
plot_df <- expand_grid(
  Phylogroup = phylo_totals$Phylogroup,
  completeness_cutoff = comp_cutoffs,
  contamination_cutoff = cont_cutoffs
) %>%
  left_join(
    phylo_counts %>%
      select(Phylogroup, completeness_cutoff, contamination_cutoff, retained),
    by = c("Phylogroup", "completeness_cutoff", "contamination_cutoff")
  ) %>%
  left_join(phylo_totals, by = "Phylogroup") %>%
  mutate(
    retained = replace_na(retained, 0),
    retained_percent = round(100 * retained / Total, 1),
    Phylogroup_label = paste0(Phylogroup, " (n=", Total, ")")
  )

p_abram_complete_genomes_phylogroup_retention_heatmap <- ggplot(plot_df,
  aes(x = factor(contamination_cutoff), y = factor(completeness_cutoff), fill = retained_percent)) +
  geom_tile() +
  geom_text(aes(label = retained_percent), size = 2.5) +
  facet_wrap(~ Phylogroup_label) +
  scale_fill_viridis_c(name = "% retained") +
  labs(
    x = "CheckM contamination cutoff (%)",
    y = "CheckM completeness cutoff (%)",
    title = "Phylogroup retention percentage across CheckM thresholds",
    subtitle = "Faceted by phylogroup") +
  theme_minimal()

print(p_abram_complete_genomes_phylogroup_retention_heatmap)


# Fisher's test to evaluate whether a phylogroup is disproportionately affected
# by combinations of completeness and contamination cutoffs

library(dplyr)
library(purrr)
library(tidyr)

phylo_df <- abram_complete_genomes %>%
  filter(
    is.finite(checkm_completeness),
    is.finite(checkm_contamination),
    !is.na(Phylogroup))

comp_cutoffs <- c(95, 97, 98, 99, 99.5)
cont_cutoffs <- c(0.5, 1, 2, 3, 5)

fisher_results <- expand_grid(
  completeness_cutoff = comp_cutoffs,
  contamination_cutoff = cont_cutoffs
) %>%
  mutate(
    results = map2(completeness_cutoff, contamination_cutoff, function(comp_cutoff, cont_cutoff) {
      
      test_df <- phylo_df %>%
        mutate(retained = checkm_completeness >= comp_cutoff & checkm_contamination <= cont_cutoff)
      
      phylogroups <- unique(test_df$Phylogroup)
      
      map_dfr(phylogroups, function(pg) {
        
        tbl <- table(
          factor(test_df$Phylogroup == pg, levels = c(FALSE, TRUE)),
          factor(test_df$retained, levels = c(FALSE, TRUE))
        )
        
        # skip cases where Fisher's test is not meaningful
        if (any(dim(tbl) < 2)) {
          return(data.frame(
            phylogroup = pg,
            odds_ratio = NA_real_,
            p_value = NA_real_
          ))
        }
        
        ft <- fisher.test(tbl)
        
        data.frame(
          phylogroup = pg,
          odds_ratio = unname(ft$estimate),
          p_value = ft$p.value
        )
      }) %>%
        mutate(
          completeness_cutoff = comp_cutoff,
          contamination_cutoff = cont_cutoff,
          p_adj = p.adjust(p_value, method = "BH")
        )
    })
  ) %>%
  select(results) %>%
  unnest(results)

# write.csv(fisher_results, "abram_et_al_phylogroup_checkm_fishers_test_raw_results.csv", row.names = FALSE)

# present Fisher's results as the number of significant phylogroups at each cutoff
# level, and the list of those phylogroups

fisher_summary_detailed <- fisher_results %>%
  filter(p_adj < 0.05) %>%
  group_by(completeness_cutoff, contamination_cutoff) %>%
  summarise(
    significant_phylogroups = n(),
    phylogroups = paste(phylogroup, collapse = ", "),
    .groups = "drop"
  )

fisher_summary_detailed

write.csv(
  fisher_summary_detailed, "abram_et_al_phylogroup_checkm_fishers_test_summary.csv", row.names = FALSE)


# plot stacked relative abundance bars, with an unfiltered "All" bar in every facet
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(purrr)

comp_cutoffs <- c(95, 97, 98, 99, 99.5)
cont_cutoffs <- c(0.5, 1, 2, 3, 5)

phylogroup_levels <- sort(unique(abram_complete_genomes$Phylogroup))

# "All" bar repeated in every contamination facet as a reference point
all_df <- expand_grid(
  completeness_cutoff = "All", contamination_cutoff = cont_cutoffs) %>%
  mutate(
    completeness_cutoff = as.character(completeness_cutoff),
    contamination_cutoff = as.character(contamination_cutoff)) %>%
  tidyr::crossing(
    abram_complete_genomes %>%
      mutate(Phylogroup = factor(Phylogroup, levels = phylogroup_levels)) %>%
      count(Phylogroup, name = "n"))

# cutoff-specific bars
cutoff_df <- expand_grid(
  completeness_cutoff = comp_cutoffs,
  contamination_cutoff = cont_cutoffs) %>%
  pmap_dfr(function(completeness_cutoff, contamination_cutoff) {
    abram_complete_genomes %>%
      filter(checkm_completeness >= completeness_cutoff,
            checkm_contamination <= contamination_cutoff) %>%
      mutate(Phylogroup = factor(Phylogroup, levels = phylogroup_levels)) %>%
      count(Phylogroup, name = "n") %>%
      mutate(completeness_cutoff = as.character(completeness_cutoff),
            contamination_cutoff = as.character(contamination_cutoff))
    })

plot_df <- bind_rows(all_df, cutoff_df) %>%
  mutate(completeness_cutoff = factor(
      completeness_cutoff,
      levels = c("All", as.character(comp_cutoffs))),
    contamination_cutoff = factor(
      contamination_cutoff,
      levels = as.character(cont_cutoffs))) %>%
  group_by(completeness_cutoff, contamination_cutoff) %>%
  mutate(relative_abundance = n / sum(n)) %>%
  ungroup()

bar_labels <- plot_df %>%
  group_by(completeness_cutoff, contamination_cutoff) %>%
  summarise(total_n = sum(n), .groups = "drop")

p_phylogroup_abundance_per_cut_off <- ggplot(plot_df,
  aes(x = completeness_cutoff, y = relative_abundance, fill = Phylogroup)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ contamination_cutoff) +
  scale_fill_brewer(palette = "Set3", drop = FALSE) +
  geom_text(data = bar_labels,
    aes(x = completeness_cutoff, y = 1.02,
      label = paste0("n = ", scales::comma(total_n))),
    inherit.aes = FALSE, size = 3, vjust = -0.5) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "CheckM completeness cutoff (%)",
    y = "Relative abundance of retained genomes",
    fill = "Phylogroup",
    title = "Phylogroup composition of the reference database across CheckM completeness/contamination cutoffs",
    subtitle = "Faceted by checkM contamination (%); 'All' bar is unfiltered and repeated in each facet for reference") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank())

print(p_phylogroup_abundance_per_cut_off)

# from assessing the above  checkM cutoffs of >= 97% completeness and <= 3% contamination will be applied for the final reference database
# This balances the number of genomes left (526) with the phylogroup representation (all phylogroups retained, with limited phylogroups being significant in the Fishers test results)
# Extract a final table based on these filters

abram_complete_genomes_checkM_cutoffs <- abram_complete_genomes %>%
  filter(
    checkm_completeness >= 97,
    checkm_contamination <= 3
  ) %>%
  arrange(Phylogroup, desc(checkm_completeness), desc(checkm_contamination))

# Save final table - these genomes are used for input to dRep for dereplication using the paramaters outlined in Supplementary Methods S2.
write.csv(
  abram_complete_genomes_checkM_cutoffs, "abram_complete_genomes_checkM_cutoffs.csv", row.names = FALSE)