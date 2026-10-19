# ============================================================
# Fig. 8. Predicted E. coli relative abundance by phylogroup
# across 42 real-world metagenomes
#
# 10 true spiked strains removed from the reference database
#
# LOW-ABUNDANCE FILTER:
# Raw strain-level results below 0.005 (0.5%) removed
# before re-normalisation.
#
# After filtering, predicted abundances are re-normalised
# WITHIN EACH TOOL AND SAMPLE so that the remaining calls
# sum to 1.
# ============================================================

# Load the necessary libraries (ensure they are installed on youe system)
library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(readxl)
library(cowplot)

# ------------------------------------------------------------
# 1. Load master dataset
# ------------------------------------------------------------
master <- read.csv("real_world_spike_ins_benchmarking_no_ground_truth_results.csv") # <-- Available on FigShare.

# ------------------------------------------------------------
# 2. Prepare abundance column prior to filtering
# ------------------------------------------------------------
master <- master %>%
  mutate(predicted_abundance = as.numeric(predicted_abundance))

# ------------------------------------------------------------
# 3. Remove low-abundance raw results
# Everything below 0.005 (0.5%) is removed.
# ------------------------------------------------------------
master <- master %>%
  filter(
    !is.na(predicted_abundance),
    predicted_abundance >= 0.005)

# ------------------------------------------------------------
# 4. Renormalise EVERY tool after filtering
# Each tool/sample combination is independently scaled
# so that the remaining strain abundances sum to 1.
#
# Note that this is a post-hoc normalisation for comparative
# visualisation, not a correction to any algorithm.
# ------------------------------------------------------------
master <- master %>%
  group_by(
    strain_profiling_tool,
    sample_id
  ) %>%
  mutate(
    predicted_abundance = predicted_abundance /
      sum(predicted_abundance, na.rm = TRUE)
  ) %>%
  ungroup()

write.csv(master, "real_world_spike_ins_no_ground_truth_combined_results_filtered_005.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. Check that renormalisation worked
# Every tool/sample combination with retained results
# should now sum to 1.
# Also identify any groups that do NOT sum to 1 (if any).
# ------------------------------------------------------------
normalisation_check <- master %>%
  group_by(
    strain_profiling_tool,
    sample_id
  ) %>%
  summarise(
    total = sum(
      predicted_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

print(normalisation_check, n = 252)

normalisation_check %>%
  filter(abs(total - 1) > 1e-6)

# ------------------------------------------------------------
# 6. Prepare main data source for figures
# ------------------------------------------------------------
master2 <- master %>%
  mutate(predicted_abundance = as.numeric(predicted_abundance),
    detected = as.logical(detected),
    phylogroup = as.character(phylogroup))

# ------------------------------------------------------------
# 7. Check for unexpected/missing phylogroups
# ------------------------------------------------------------
print(
  table(
    master2$phylogroup,
    useNA = "ifany"))

# ------------------------------------------------------------
# 8. Calculate predicted abundance at phylogroup level
# ------------------------------------------------------------
phylogroup_abundance <- master2 %>%
  group_by(
    strain_profiling_tool,
    sample_id,
    phylogroup
  ) %>%
  summarise(
    predicted_abundance = sum(
      predicted_abundance,
      na.rm = TRUE
    ),
    .groups = "drop")

# ------------------------------------------------------------
# 9. Count detected phylogroups and strains per sample/tool
# ------------------------------------------------------------
detected_counts <- master2 %>%
  filter(
    detected == TRUE
  ) %>%
  group_by(
    strain_profiling_tool,
    sample_id
  ) %>%
  summarise(
    n_phylogroups = n_distinct(phylogroup),
    n_strains = n_distinct(strain_id),
    .groups = "drop"
  ) %>%
  mutate(
    detection_label = paste0(
      n_phylogroups,
      "(",
      n_strains,
      ")"
    )
  )

# ------------------------------------------------------------
# 10. Define plotting data
# ------------------------------------------------------------
plot_data <- phylogroup_abundance

# ------------------------------------------------------------
# 11. Define phylogroup order and colours
# ------------------------------------------------------------
phylo_levels <- c(
  "A",
  "B1",
  "B2-1",
  "B2-2",
  "C",
  "D1",
  "D2",
  "D3",
  "E1",
  "E2",
  "F",
  "G")

phylo_colours <- setNames(
  RColorBrewer::brewer.pal(
    12,
    "Paired"
  ),
  phylo_levels)

plot_data <- plot_data %>%
  mutate(
    phylogroup = factor(
      phylogroup,
      levels = phylo_levels))

# ------------------------------------------------------------
# 12. Define sample order
# ------------------------------------------------------------
sample_levels <- master2 %>%
  distinct(sample_id) %>%
  arrange(sample_id) %>%
  pull(sample_id)

plot_data <- plot_data %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_levels))

detected_counts <- detected_counts %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_levels))

# ------------------------------------------------------------
# 13. Define tool order
# ------------------------------------------------------------
tool_levels <- sort(
  unique(
    as.character(
      master2$strain_profiling_tool)))

plot_data <- plot_data %>%
  mutate(
    strain_profiling_tool = factor(
      strain_profiling_tool,
      levels = tool_levels))

detected_counts <- detected_counts %>%
  mutate(
    strain_profiling_tool = factor(
      strain_profiling_tool,
      levels = tool_levels))

# ------------------------------------------------------------
# 14. Main stacked bar plot
# ------------------------------------------------------------
p_phylogroup_removed_strains <- ggplot(plot_data,
  aes(x = sample_id,y = predicted_abundance, fill = phylogroup)) +
  geom_col(
    color = "black",linewidth = 0.2, width = 0.8,
    position = position_stack(reverse = FALSE)) +
  # Number of detected reference strains above each bar
  geom_text(
    data = detected_counts,
    aes(x = sample_id, y = 1.02, label = detection_label),
    inherit.aes = FALSE, size = 2.1, vjust = -0.1) +
  # Abundance labels for larger segments
  geom_text(
    aes(label = ifelse(
        predicted_abundance >= 0.25,
        sprintf("%.2f", predicted_abundance), "")),
    position = position_stack(vjust = 0.5),
        size = 2.75, color = "white", fontface = "bold", angle = 0, na.rm = TRUE) +
  # One row per tool
  facet_grid(
    rows = vars(strain_profiling_tool),
    scales = "fixed",
    labeller = labeller(
      strain_profiling_tool = function(x) {
        stringr::str_wrap(
          x, width = 15)})) +
  scale_fill_manual(
    values = phylo_colours, breaks = phylo_levels,
    drop = FALSE, name = expression(
      italic("E. coli") ~ "phylogroup")) +
  # Keep y-axis on the 0-1 scale
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1.08),
    expand = expansion(
      mult = c(0, 0))) +
  labs(
    x = "Real-world metagenome sample",
    y = "Predicted relative abundance") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    axis.title.y = element_text(margin = margin(r = 10), size = 11),
    strip.position = "right",
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 11),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10,face = "bold"))

print(p_phylogroup_removed_strains)

# ------------------------------------------------------------
# 15. Extract legend for custom figure layout
# ------------------------------------------------------------
phylogroup_legend <- get_legend(p_phylogroup_removed_strains)

# ------------------------------------------------------------
# 16. Save main filtered figure
# ------------------------------------------------------------
ggsave("real_world_spike_ins_no_ground_truth_phylogroup_level_stacked_bars_abundance_filtered_005.pdf",
  plot = p_phylogroup_removed_strains, width = 18, height = 14)

# ============================================================
# 17. Ground-truth bar (added to figure for reference)
# ============================================================
ground_truth_phylo <- data.frame(
  phylogroup = phylo_levels,
  
  abundance = c(
    0.12,  # A
    0.08,  # B1
    0.20,  # B2-1
    0.50,  # B2-2
    0.00,  # C
    0.00,  # D1
    0.00,  # D2
    0.00,  # D3
    0.00,  # E1
    0.06,  # E2
    0.04,  # F
    0.00   # G
  )
)

# ------------------------------------------------------------
# 18. Keep only phylogroups represented in ground truth
# ------------------------------------------------------------
ground_truth_phylo <- ground_truth_phylo %>%
  filter(
    abundance > 0
  ) %>%
  mutate(
    phylogroup = factor(
      phylogroup,
      levels = phylo_levels))

# ------------------------------------------------------------
# 19. Calculate label positions to match stacking order
# ------------------------------------------------------------
ground_truth_phylo <- ground_truth_phylo %>%
  arrange(
    desc(phylogroup)
  ) %>%
  mutate(
    label_y = cumsum(
      abundance
    ) - abundance / 2,
    ground_truth_label = paste0(
      as.character(phylogroup),
      ", ",
      sprintf(
        "%.2f",
        abundance
      )
    ),
    label_colour = ifelse(
      phylogroup %in% c("F"),
      "black",
      "white"))

# ------------------------------------------------------------
# 20. Create vertical stacked ground-truth bar
# ------------------------------------------------------------
ground_truth_bar <- ggplot(ground_truth_phylo,
  aes(x = "Ground truth", y = abundance, fill = phylogroup)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  geom_text(
    aes(x = "Ground truth", y = label_y,
      label = ground_truth_label, color = label_colour),
    size = 2.5, fontface = "bold", show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_manual(
    values = phylo_colours, breaks = phylo_levels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))) +
  theme_void() +
  theme(legend.position = "none")

# ------------------------------------------------------------
# 21. Add label beneath ground-truth bar
# ------------------------------------------------------------
ground_truth_label <- ggdraw() +
  draw_label(
    "Ground truth\n(expected\nphylogroup\nproportions)",
    fontface = "bold", size = 10, x = 0.5, y = 0.35)

ground_truth_panel <- plot_grid(
  ground_truth_bar,
  ground_truth_label,
  ncol = 1,
  rel_heights = c(1, 0.25))

# ------------------------------------------------------------
# 22. Position ground-truth bar underneath legend
# ------------------------------------------------------------
right_hand_panel <- ggdraw() +
  # Existing phylogroup legend
  draw_plot(
    phylogroup_legend,
    x = 0, y = 0.43,
    width = 1, height = 0.57) +
  # Ground-truth bar
  draw_plot(
    ground_truth_panel,
    x = 0.275, y = 0.22,
    width = 0.35, height = 0.32)

# ------------------------------------------------------------
# 23. Combine main plot with custom legend area
# ------------------------------------------------------------
p_phylogroup_removed_strains_with_truth <- plot_grid(
  p_phylogroup_removed_strains +
    theme(legend.position = "none"),
  right_hand_panel,
  ncol = 2,
  rel_widths = c(1, 0.18))

print(p_phylogroup_removed_strains_with_truth)

# ------------------------------------------------------------
# 24. Save final figure
# ------------------------------------------------------------
ggsave("Fig8_abundance_filtered_005.pdf", plot = p_phylogroup_removed_strains_with_truth, width = 18, height = 14)


# ============================================================
# Table S9 - extraction and plotting script
# Genomic relatedness of PathoScope and StrainGE predictions
# to removed ground-truth strains
#
# For each removed ground-truth strain and each tool:
#   1. Find detected reference strains belonging to the same
#      phylogroup as the removed truth
#   2. Identify the detected reference with the highest pairwise ANI
#   3. Count the number of samples in which that strain was detected
#   4. Record its ANI rank among all 294 available references
# ============================================================

# (Re)load libraries
library(dplyr)
library(readr)
library(stringr)

# ============================================================
# 1. Files (all available - see Data summary)
# ============================================================
ANI_FILE <- "fastANI_dRep_genomes_all_vs_all_real_world_spike_ins.tsv" # <-- Available on FigShare.
METADATA_FILE <- "304_ecoli_strains_reference_genomes.csv" # <-- Available on FigShare.
MASTER_FILE <- "real_world_spike_ins_no_ground_truth_combined_results.csv" # <-- Available on FigShare.
OUTPUT_FILE <- "Table_S9_PathoScope_StrGE_removed_truth_ANI_summary.tsv" # <-- Available on FigShare.

# ============================================================
# 2. Helper function
# Extract base GCA accession to match output of tools, e.g.,
# GCA_000493595.1_EcoPMV1_genomic.fna -> GCA_000493595
# ============================================================
get_gca <- function(x) {
  x <- basename(as.character(x))
  out <- str_extract(
    x,
    "GCA_[0-9]+"
  )
  ifelse(
    is.na(out),
    NA_character_,
    out
  )
}

# ============================================================
# 3. Read master benchmarking results
# ============================================================
master <- read_csv(MASTER_FILE,
                   show_col_types = FALSE)
# Make sure important columns have correct types
master <- master %>%
  mutate(
    strain_id = as.character(strain_id),
    sample_id = as.character(sample_id),
    strain_profiling_tool = as.character(strain_profiling_tool),
    detected = as.logical(detected))

# ============================================================
# 4. Read fastANI results
# ============================================================
ani <- read_tsv(ANI_FILE,
                col_names = FALSE, show_col_types = FALSE)

colnames(ani) <- c(
  "query",
  "reference",
  "ANI",
  "matching_fragments",
  "total_fragments")

ani <- ani %>%
  mutate(
    query_gca = get_gca(query),
    reference_gca = get_gca(reference),
    ANI = as.numeric(ANI))

# ============================================================
# 5. Define the 10 removed ground-truth strains
# ============================================================
true_strains <- c(
  "GCA_000493595.1_EcoPMV1_genomic",
  "GCA_000025745.1_ASM2574v1_genomic",
  "GCA_001280345.1_ASM128034v1_genomic",
  "GCA_001900395.1_ASM190039v1_genomic",
  "GCA_002863845.1_ASM286384v1_genomic",
  "GCA_002996945.1_ASM299694v1_genomic",
  "GCA_002812545.1_ASM281254v1_genomic",
  "GCA_003073815.1_ASM307381v1_genomic",
  "GCA_000245515.1_ASM24551v1_genomic",
  "GCA_002948655.1_ASM294865v1_genomic"
) %>%
  get_gca()

# ============================================================
# 6. Read phylogroup metadata
# ============================================================
metadata <- read_csv(METADATA_FILE,
                     show_col_types = FALSE)

metadata <- metadata %>%
  transmute(
    candidate_gca = get_gca(accession),
    phylogroup = as.character(Phylogroup)
  ) %>%
  distinct(candidate_gca, .keep_all = TRUE)

# ============================================================
# 7. Define the 294 available reference genomes
# All genomes represented in fastANI minus  10 removed true strains
# ============================================================
matrix_genomes <- union(
  ani$query_gca,
  ani$reference_gca)

candidate_genomes <- setdiff(
  matrix_genomes,
  true_strains)

cat(
  "Available candidate genomes:",
  length(candidate_genomes),
  "\n")

# ============================================================
# 8. Generate complete ANI rankings
# One row = one removed truth × one available reference
# ============================================================
all_rankings <- lapply(
  true_strains,
  function(true_strain) {
    
    # --------------------------------------------------------
    # ANI where true strain is query
    # --------------------------------------------------------
    q <- ani %>%
      filter(query_gca == true_strain) %>%
      transmute(
        candidate_gca = reference_gca,
        ANI,
        matching_fragments,
        total_fragments)
    
    # --------------------------------------------------------
    # ANI where true strain is reference
    # --------------------------------------------------------
    r <- ani %>%
      filter(reference_gca == true_strain) %>%
      transmute(
        candidate_gca = query_gca,
        ANI,
        matching_fragments,
        total_fragments)
    
    # --------------------------------------------------------
    # Combine both orientations
    # --------------------------------------------------------
    tmp <- bind_rows(q, r) %>%
      filter(
        candidate_gca != true_strain,
        !candidate_gca %in% true_strains
      ) %>%
      arrange(desc(ANI)) %>%
      distinct(
        candidate_gca,
        .keep_all = TRUE
      ) %>%
      arrange(desc(ANI))
    
    # --------------------------------------------------------
    # Global ANI rank among available 294 references
    # --------------------------------------------------------
    tmp <- tmp %>%
      mutate(
        ANI_rank = row_number(),
        true_strain = true_strain
      ) %>%
      left_join(
        metadata,
        by = "candidate_gca")
    
    # --------------------------------------------------------
    # Add true-strain phylogroup
    # --------------------------------------------------------
    true_phylogroup <- metadata %>%
      filter(candidate_gca == true_strain) %>%
      pull(phylogroup)
    
    tmp %>%
      mutate(
        true_phylogroup = true_phylogroup)
  }
) %>%
  bind_rows()

# ============================================================
# 9. Identify unique detected strains for each tool
# ============================================================
detected_strains <- master %>%
  filter(
    strain_profiling_tool %in% c(
      "PathoScope",
      "StrainGE"
    ),
    !is.na(predicted_abundance),
    predicted_abundance >= 0.005
  ) %>%
  distinct(
    strain_profiling_tool,
    strain_id
  ) %>%
  mutate(
    candidate_gca = get_gca(strain_id))

# ============================================================
# 10. Add detection status for each tool
# ============================================================
all_rankings <- all_rankings %>%
  mutate(
    detected_by_tool = FALSE)

# ============================================================
# 11. Function to generate one Table S9 row
# for a given tool + removed truth
# ============================================================
make_tool_summary <- function(tool_name) {
  
  # All strains detected by this tool
  tool_strains <- detected_strains %>%
    filter(
      strain_profiling_tool == tool_name
    ) %>%
    pull(candidate_gca)
  # ----------------------------------------------------------
  # For each removed truth:
  # find detected candidates in the same phylogroup
  # ----------------------------------------------------------
  summaries <- lapply(
    true_strains,
    function(true_strain) {
      
      true_info <- all_rankings %>%
        filter(
          true_strain == !!true_strain
        ) %>%
        slice(1) %>%
        select(
          true_strain,
          true_phylogroup
        )
      
      true_phylogroup <- true_info$true_phylogroup
      
      # ------------------------------------------------------
      # Candidate predictions from this tool
      # within the same phylogroup
      # ------------------------------------------------------
      candidates <- all_rankings %>%
        filter(
          true_strain == !!true_strain,
          true_phylogroup == !!true_phylogroup,
          phylogroup == !!true_phylogroup,
          candidate_gca %in% tool_strains
        ) %>%
        arrange(desc(ANI))
      
      # ------------------------------------------------------
      # No candidate in same phylogroup
      # ------------------------------------------------------
      if (nrow(candidates) == 0) {
        
        return(
          tibble(
            Tool = tool_name,
            `Removed ground-truth strain` = true_strain,
            Phylogroup = true_phylogroup,
            `Closest predicted strain` = NA_character_,
            `Samples detected in (n = 42)` = 0,
            `Predicted phylogroup` = NA_character_,
            `Pairwise ANI (%)` = NA_real_,
            `ANI rank among 294 references` = NA_real_
          )
        )
      }
      
      # ------------------------------------------------------
      # Closest detected reference
      # ------------------------------------------------------
      best <- candidates %>%
        slice(1)
      
      # ------------------------------------------------------
      # Identify actual strain name / ID
      #
      # We try to use strain_id from master.
      # If multiple rows exist, take the first.
      # ------------------------------------------------------
      display_name <- master %>%
        filter(
          strain_profiling_tool == tool_name,
          detected == TRUE,
          get_gca(strain_id) == best$candidate_gca
        ) %>%
        pull(strain_id) %>%
        unique()
      
      display_name <- if (
        length(display_name) > 0
      ) {
        display_name[1]
      } else {
        best$candidate_gca
      }
      
      # ------------------------------------------------------
      # Number of samples in which this strain was detected
      # ------------------------------------------------------
      n_samples <- master %>%
        filter(
          strain_profiling_tool == tool_name,
          get_gca(strain_id) == best$candidate_gca,
          !is.na(predicted_abundance),
          predicted_abundance >= 0.005
        ) %>%
        summarise(
          n = n_distinct(sample_id)
        ) %>%
        pull(n)
      
      # ------------------------------------------------------
      # Return final row
      # ------------------------------------------------------
      tibble(
        Tool = tool_name,
        `Removed ground-truth strain` = true_strain,
        Phylogroup = true_phylogroup,
        `Closest predicted strain` = display_name,
        `Samples detected in (n = 42)` = n_samples,
        `Predicted phylogroup` = best$phylogroup,
        `Pairwise ANI (%)` = round(best$ANI, 3),
        `ANI rank among 294 references` = best$ANI_rank
      )
    }
  )
  
  bind_rows(summaries)
}


# ============================================================
# 12. Generate summaries for both tools
# ============================================================
pathoscope_summary <- make_tool_summary("PathoScope")
strainge_summary <- make_tool_summary("StrainGE")

# ============================================================
# 13. Combine both tools
# ============================================================
table_s9 <- bind_rows(
  strainge_summary,
  pathoscope_summary)

# ============================================================
# 14. Order by phylogroup
# ============================================================
phylogroup_order <- c(
  "A",
  "B1",
  "B2-1",
  "B2-2",
  "E2",
  "F")

table_s9 <- table_s9 %>%
  mutate(
    Phylogroup = factor(
      Phylogroup,
      levels = phylogroup_order
    )
  ) %>%
  arrange(
    Phylogroup,
    `Removed ground-truth strain`,
    Tool
  ) %>%
  mutate(
    Phylogroup = as.character(Phylogroup))

# ============================================================
# 15. Print and save
# ============================================================
print(table_s9)
write_tsv(table_s9, OUTPUT_FILE)
cat("\nTable S9 saved to:", OUTPUT_FILE,"\n")