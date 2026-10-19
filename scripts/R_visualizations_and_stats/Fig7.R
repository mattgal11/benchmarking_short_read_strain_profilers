# ============================================================
# Fig. 7: Predicted E. coli strain relative abundance across 
# 42 real-world metagenomes with identical 10-strain spike-in
# ============================================================

# ============================================================
# Fig. 7a: Predicted E. coli strain relative abundance across
# 42 real-world metagenomes with identical 10-strain spike-in
# One row = strain-profiling tool
# One bar = metagenome sample
# Stack = 10 true strains + "Other"
# True abundance is shown in the legend
# ============================================================

# ------------------------------------------------------------
# 1. Load libraries (ensure they are installed on your system)
# ------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(readxl)
library(grid)
library(ggtext)
library(cowplot)

# ------------------------------------------------------------
# 2. Load master final benchmarking dataset
# ------------------------------------------------------------
master <- read.csv("real_world_spike_ins_benchmarking_results.csv") # <-- Available on FigShare

# ------------------------------------------------------------
# 2b. Normalise Strainify abundance estimates to sum to 1
# Strainify reports abundances on a 0-100 scale and, for these
# benchmark samples, the raw estimates do not sum to 100 due 
# to a separatley reported solver issue. Note: this is post-
# hoc normalisation for comparative visualisation, not a 
# correction to the Strainify algorithm
# ------------------------------------------------------------
master <- master %>%
  group_by(strain_profiling_tool, sample_id) %>%
  mutate(
    predicted_abundance = if_else(
      strain_profiling_tool == "Strainify",
      predicted_abundance / sum(predicted_abundance, na.rm = TRUE),
      predicted_abundance
    )
  ) %>%
  ungroup()

# check the normalisation worked
master %>%
  filter(strain_profiling_tool == "Strainify") %>%
  group_by(sample_id) %>%
  summarise(
    total = sum(predicted_abundance, na.rm = TRUE))

# ------------------------------------------------------------
# 3. Define the 10 true spiked strains
# ------------------------------------------------------------
# IDs correspond to the strain_id column in master dataset
spiked_strains <- c(
  "GCA_000493595" = "PMV-1",
  "GCA_000025745" = "IHE3034",
  "GCA_001280345" = "SF-468",
  "GCA_001900395" = "D8",
  "GCA_002863845" = "STEC299",
  "GCA_002996945" = "AR_0081",
  "GCA_002812545" = "SC516",
  "GCA_003073815" = "AR437",
  "GCA_000245515" = "RM12579",
  "GCA_002948655" = "AR_0077")

# ------------------------------------------------------------
# 4. Define true relative abundances of spiked strains
# ------------------------------------------------------------
true_abundances <- c(
  "PMV-1"   = 0.10,
  "IHE3034" = 0.10,
  "SF-468"  = 0.10,
  "D8"      = 0.10,
  "STEC299" = 0.10,
  "AR_0081" = 0.20,
  "SC516"   = 0.12,
  "AR437"   = 0.08,
  "RM12579" = 0.06,
  "AR_0077" = 0.04)

# The first 5 strains are all within phylogroup B2-2, are increasignly genetically related (%ANI), and present in equal abundance
# The last 5 strains have more diverse phylogroups and have differential abundances

# ------------------------------------------------------------
# 5. Prepare master dataset
# ------------------------------------------------------------
master2 <- master %>%
  mutate(
    predicted_abundance = as.numeric(predicted_abundance),
    
    # Assign each reference genome to one of the 10 true strains.
    # All other E. coli genomes are grouped into "Other".
    strain_category = case_when(
      strain_id %in% names(spiked_strains) ~
        unname(spiked_strains[strain_id]),
      TRUE ~ "Other"
    )
  )

# ------------------------------------------------------------
# 6. Collapse predicted abundance into strain categories
# ------------------------------------------------------------
# For each tool and metagenome, sum the predicted abundance assigned to each true strain.
# Predictions assigned to any other E. coli reference genome are collapsed into the "Other" category.

plot_data <- master2 %>%
  group_by(
    strain_profiling_tool,
    sample_id,
    strain_category
  ) %>%
  summarise(
    predicted_relative_abundance = if (
      all(is.na(predicted_abundance))
    ) {
      NA_real_
    } else {
      sum(predicted_abundance, na.rm = TRUE)
    },
    .groups = "drop"
  )

# ------------------------------------------------------------
# 7. Define sample order
# ------------------------------------------------------------
# Alphabetical sample order ensures every tool has the same ordering of the 42 metagenomes.
# In this case, all ERSxxx followed by all SRSxxx metagenomes

sample_levels <- master2 %>%
  distinct(sample_id) %>%
  arrange(sample_id) %>%
  pull(sample_id)

# ------------------------------------------------------------
# 8. Define tool order
# ------------------------------------------------------------
tool_levels <- sort(
  unique(as.character(master2$strain_profiling_tool)))

# ------------------------------------------------------------
# 9. Define strain stacking / legend order
# ------------------------------------------------------------
# First five strains represent the equal-abundance B2-2 strains.
# Following five represent the differential-abundance / phylogroup strains.
# "Other" represents abundance assigned to non-spiked E. coli reference genomes present in the 304 E. coli reference genome database

strain_levels <- c(
  "Other",
  "PMV-1",
  "IHE3034",
  "SF-468",
  "D8",
  "STEC299",
  "AR_0081",
  "SC516",
  "AR437",
  "RM12579",
  "AR_0077")

# ------------------------------------------------------------
# 10. Apply factor levels
# ------------------------------------------------------------
plot_data <- plot_data %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_levels
    ),
    strain_profiling_tool = factor(
      strain_profiling_tool,
      levels = tool_levels
    ),
    strain_category = factor(
      strain_category,
      levels = strain_levels))

# ------------------------------------------------------------
# 11. Plot, print and save
# ------------------------------------------------------------
p_final_benchmark_stacked <- ggplot(plot_data,
    aes(x = sample_id, y = predicted_relative_abundance, fill = strain_category)) +
  geom_col(color = "black", linewidth = 0.2, width = 0.85,
    position = position_stack(reverse = FALSE), na.rm = TRUE) +
  geom_text(aes(label = ifelse(
        predicted_relative_abundance >= 0.1,
        sprintf("%.2f", predicted_relative_abundance), "")),
       position = position_stack(vjust = 0.5),
       size = 2.5, color = "white", fontface = "bold", angle = 0, na.rm = TRUE) +
  facet_grid(rows = vars(strain_profiling_tool),
    scales = "fixed",
    labeller = labeller(
      strain_profiling_tool = function(x) {
        stringr::str_wrap(x, width = 15)
      })) +
  scale_fill_brewer(
    palette = "Paired",
    name = "<b><i>E. coli</i> strain</b><br><b>(phylogroup;</b><br><b>true abundance)</b>",
    breaks = strain_levels,
    labels = c(
      "Other"   = "Other",
      "PMV-1"   = "PMV-1 (B2-2; 0.10)",
      "IHE3034" = "IHE3034 (B2-2; 0.10)",
      "SF-468"  = "SF-468 (B2-2; 0.10)",
      "D8"      = "D8 (B2-2; 0.10)",
      "STEC299" = "STEC299 (B2-2; 0.10)",
      "AR_0081" = "AR_0081 (B2-1; 0.20)",
      "SC516"   = "SC516 (A; 0.12)",
      "AR437"   = "AR437 (B1; 0.08)",
      "RM12579" = "RM12579 (E2; 0.06)",
      "AR_0077" = "AR_0077 (F; 0.04)"
    )) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.01))) +
  labs(
    x = "Real-world spiked metagenome sample",
    y = "Predicted relative abundance") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 6),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text.x.bottom = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold", size = 10, angle = 0),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "right",
    legend.text = element_text(size = 6.75),
    legend.title = ggtext::element_markdown(size = 8, face = "bold"))

print(p_final_benchmark_stacked)
# ggsave("Fig7a.pdf", plot = p_final_benchmark_stacked, width = 18, height = 12)

# ------------------------------------------------------------
# 12. Add representative ground-truth bar to legend
# Load cowplot library and create ground truth data
# ------------------------------------------------------------
ground_truth_data <- data.frame(
  strain_category = names(true_abundances),
  abundance = as.numeric(true_abundances))

# Add Other = 0 so categories match the main plot
ground_truth_data <- bind_rows(
  data.frame(
    strain_category = "Other",
    abundance = 0),
  ground_truth_data)

# Use the same strain order as the main plot and
# calculate the midpoint of each segment for labels
ground_truth_data <- ground_truth_data %>%
  mutate(
    strain_category = factor(
      strain_category,
      levels = strain_levels),
    label_y = rev(cumsum(rev(abundance)) - rev(abundance) / 2),
    ground_truth_label = ifelse(
      abundance > 0,
      paste0(
        as.character(strain_category),
        ", ",
        sprintf("%.2f", abundance)),
      ""),
    # Use black text for AR_0077 because of its light fill colour
    label_colour = ifelse(
      strain_category == "AR_0077",
      "black",
      "white"))

# ------------------------------------------------------------
# 13. Create one vertical stacked ground-truth bar
# ------------------------------------------------------------
ground_truth_bar <- ggplot(ground_truth_data,
                           aes(x = "Ground truth", y = abundance, fill = strain_category)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.2,
    position = position_stack(reverse = FALSE)) +
   geom_text(aes(y = label_y,
      label = ground_truth_label,
      color = label_colour),
    size = 2.25, fontface = "bold",
    lineheight = 0.9, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_brewer(
    palette = "Paired",
    breaks = strain_levels, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0))) +
  theme_void() +
  theme(legend.position = "none")

# ------------------------------------------------------------
# 14. Extract existing strain legend from main plot
# ------------------------------------------------------------
strain_legend <- cowplot::get_legend(
  p_final_benchmark_stacked +
    theme(legend.position = "right"))

# ------------------------------------------------------------
# 15. Add "Ground truth" label beneath bar
# ------------------------------------------------------------
ground_truth_label <- ggdraw() +
  draw_label(
    "Ground truth \n(expected \nproportions)",
    fontface = "bold",
    size = 8, x = 0.5, y = 0.25)

# Combine ground-truth bar and label
ground_truth_panel <- plot_grid(
  ground_truth_bar,
  ground_truth_label,
  ncol = 1,
  rel_heights = c(1, 0.15))

# ------------------------------------------------------------
# 16. Position and combine  existing legend + ground-truth bar
# ------------------------------------------------------------
right_hand_panel <- ggdraw() +
  # Existing strain legend
  draw_plot(
    strain_legend,
    x = 0, y = 0.45,
    width = 1, height = 0.55) +
  # Ground-truth bar
  draw_plot(
    ground_truth_panel,
    x = 0.1, y = 0.225,
    width = 0.7, height = 0.30)

# ------------------------------------------------------------
# 17. Add custom legend area to main plot
# ------------------------------------------------------------
p_final_benchmark_stacked_with_truth <- plot_grid(
  p_final_benchmark_stacked +
    theme(legend.position = "none"),
  right_hand_panel,
  ncol = 2,
  rel_widths = c(1, 0.18))

print(p_final_benchmark_stacked_with_truth)

# ------------------------------------------------------------
# 18. Save final figure
# ------------------------------------------------------------
ggsave("Fig7a_with_ground_truth.pdf", 
       plot = p_final_benchmark_stacked_with_truth, width = 18, height = 12)


##############################################################

#=============================================================
# Fig. 7b: Mean absolute proportional error per tool
# Second part of this plot, collapsing tools to mean APE +/ SD
#=============================================================

# ------------------------------------------------------------
# 1. Restrict to the 10 true spiked strains with valid APE
# ------------------------------------------------------------
df_spiked_APE <- master2 %>%
  filter(
    strain_id %in% names(spiked_strains),
    strain_present_in_sample == TRUE,
    detected == TRUE
  ) %>%
  mutate(
    absolute_proportional_error =
      as.numeric(absolute_proportional_error)
  ) %>%
  filter(
    !is.na(absolute_proportional_error))

# ------------------------------------------------------------
# 2. Calculate mean APE across true strains within each sample
# ------------------------------------------------------------
df_sample_mean_APE <- df_spiked_APE %>%
  group_by(
    sample_id,
    strain_profiling_tool
  ) %>%
  summarise(
    mean_APE_sample = mean(
      absolute_proportional_error,
      na.rm = TRUE
    ),
    .groups = "drop")

# ------------------------------------------------------------
# 3. Calculate mean +/- SD across the 42 metagenomes per tool
# ------------------------------------------------------------
df_profiler_summary_APE <- df_sample_mean_APE %>%
  group_by(strain_profiling_tool) %>%
  summarise(
    mean_APE = mean(
      mean_APE_sample,
      na.rm = TRUE
    ),
    sd_APE = sd(
      mean_APE_sample,
      na.rm = TRUE
    ),
    n_samples = n(),
    .groups = "drop")

print(df_profiler_summary_APE)

# ------------------------------------------------------------
# 4. Manually set tool colours and then plot, print and save
# ------------------------------------------------------------
tool_colours <- c(
  "PanTax"                       = "#1B9E77",
  "PathoScope"                   = "#D95F02",
  "StrainGE"                     = "#7570B3",
  "Strainify"                    = "#E7298A",
  "StrainR2"                     = "#66A61E",
  "StrainScan"                   = "#E6AB02",
  "StrainScan (super low depth)" = "#A6761D"
  )

p_final_benchmark_mean_APE_bar <- ggplot(df_profiler_summary_APE,
    aes(x = strain_profiling_tool,
        y = mean_APE)) +
  geom_col(aes(fill = strain_profiling_tool),
    width = 0.7) +
  geom_errorbar(
    aes(ymin = pmax(0, mean_APE - sd_APE),
        ymax = mean_APE + sd_APE), width = 0.2) +
  geom_text(
    aes(y = mean_APE + sd_APE,
        label = sprintf("%.2f", mean_APE)),
        vjust = -0.6, size = 5) +
  coord_cartesian(ylim = c(0, 3.88)) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 15)}) +
  scale_y_continuous(
    name = "Mean absolute proportional error",
    expand = expansion(mult = c(0, 0.0))) +
  scale_fill_manual(values = tool_colours) +
  labs(x = "Strain-level profiling tool") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
        axis.title.x = element_text(vjust = -0.75))

print(p_final_benchmark_mean_APE_bar)
ggsave("Fig7b.pdf", plot = p_final_benchmark_mean_APE_bar, width = 12, height = 6)


# ============================================================
# Cowplot combined Fig. 7 for final publication
# ============================================================

p7a <- p_final_benchmark_stacked_with_truth
p7b <- p_final_benchmark_mean_APE_bar

plot_grid(p7a, p7b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(2.75, 1))
ggsave("combined_Fig7.pdf", width = 16, height = 10)