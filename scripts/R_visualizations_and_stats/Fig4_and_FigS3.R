# Final plotting scripts for Fig. 4 and Fig. S3 (part of the same results narrative)

### Fig. 4a - Predicted relative abundance across simulated metagenome samples ###

# Load the required libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)

# ---- user-specific strain IDs (from NCBI reference download)
k12_id   <- "GCF_000005845.2_ASM584v2_genomic"
sakai_id <- "GCF_000008865.2_ASM886v2_genomic"

# ---- helper build labels
build_proportion_label <- function(sid) {
  sid <- str_remove(sid, "^HiSeq_")
  sid <- str_remove(sid, "_rep[0-9]+$")
  sid <- str_remove(sid, "^(20M|50M|100M)_")
  sid <- str_remove_all(sid, "K12|Sakai|k12|sakai|K12_|Sakai_")
  
  nums <- str_extract_all(sid, "\\d*\\.?\\d+")[[1]]
  if (length(nums) < 2) return(sid)
  
  p1 <- as.numeric(nums[length(nums) - 1])
  p2 <- as.numeric(nums[length(nums)])
  
  if (!is.na(p1) && p1 > 1) p1 <- p1 / 100
  if (!is.na(p2) && p2 > 1) p2 <- p2 / 100
  
  sprintf("%.2f:%.2f", p1, p2)
}

# ---- 1) Prepare master copy
master <- read.csv("Simulated_metagenomes_all_references_raw_predicted_abundances.csv") # Available on FigShare
master2 <- master %>%
  mutate(
    sample_clean = sub("^HiSeq_", "", sample_id),
    sample_base  = sub("_rep[0-9]+$", "", sample_clean),
    depth        = ifelse(!is.na(depth_of_sequencing) & depth_of_sequencing != "",
                          as.character(depth_of_sequencing),
                          stringr::str_extract(sample_base, "20M|50M|100M")),
    sample_label = vapply(as.character(sample_base), build_proportion_label, character(1)),
    sample_key   = ifelse(!is.na(depth),
                          paste0(depth, "_", sample_label),
                          paste0("noDepth_", sample_label)),
    strain_category = dplyr::case_when(
      strain_id == k12_id   ~ "K12-MG1655",
      strain_id == sakai_id ~ "O157:H7 Sakai",
      TRUE                  ~ "Other"
    )
  )

# ---- 2) Collapse predicted_abundance per replicate into strain categories (sum per replicate)
per_rep_category <- master2 %>%
  group_by(strain_profiling_tool, sample_with_rep = sample_clean, sample_key, sample_label, strain_category) %>%
  summarise(
    sum_abundance = if (all(is.na(predicted_abundance))) NA_real_ else sum(predicted_abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 3) Average across replicates (mean of replicate-level sums), keep sample_key so depth separation remains
mean_abundance <- per_rep_category %>%
  group_by(strain_profiling_tool, sample_key, sample_label, strain_category) %>%
  summarise(
    mean_predicted_abundance = mean(sum_abundance, na.rm = TRUE),
    n_reps_contributed = sum(!is.na(sum_abundance)),
    .groups = "drop"
  )

# ---- 4) Ensure stacking order: Other (top) -> K12 -> Sakai
mean_abundance <- mean_abundance %>%
  mutate(strain_category = factor(strain_category, levels = c("Other", "K12-MG1655", "O157:H7 Sakai")))

# ---- 5) Extract depth and set ordering for labels
mean_abundance <- mean_abundance %>%
  mutate(depth = str_extract(sample_key, "20M|50M|100M"))  # depth is inside sample_key (or NA)

depth_levels <- c("20M", "50M", "100M")
mean_abundance <- mean_abundance %>% mutate(depth = factor(depth, levels = depth_levels))

# Build consistent order for sample_label across the whole figure:
# order by depth then by sample_label so labels group logically
sample_order_tbl <- mean_abundance %>%
  distinct(depth, sample_label) %>%
  arrange(factor(depth, levels = depth_levels), sample_label)
sample_label_levels <- unique(sample_order_tbl$sample_label)
mean_abundance <- mean_abundance %>% mutate(sample_label = factor(sample_label, levels = sample_label_levels))

# ---- 6) PLOTTING: use sample_label as the x aesthetic (shared across depths)
# Use scales = "fixed" so x categories are identical across depth rows -> aligned columns
p_stacked_depth <- ggplot(mean_abundance,
                          aes(x = sample_label,
                              y = mean_predicted_abundance,
                              fill = strain_category)) +
  geom_col(color = "black", width = 0.8, position = position_stack(reverse = FALSE)) +
  geom_text(aes(label = ifelse(mean_predicted_abundance > 0.11, sprintf("%.2f", mean_predicted_abundance), "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, color = "white", fontface = "bold", angle = 90) +
  facet_grid(rows = vars(depth), cols = vars(strain_profiling_tool), 
             labeller = labeller(
               strain_profiling_tool = function(x)
                 stringr::str_wrap(x, width = 16)),scales = "fixed") +
  scale_fill_brewer(palette = "Dark2",
                    name = expression(italic("E. coli") ~ "strain")) +
  labs(
    #title = "Predicted relative abundance per sample per strain-level , faceted by sequencing effort",
    #subtitle = "'Other' includes all additional strains within reference genome database; all values are means across triplicates; no threshold applied.",
    x = "True K12-MG1655 : O157:H7 str. Sakai ratio in simulated metagenome sample",
    y = "Predicted relative abundance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    # rotate labels; reduce size to mitigate crowding
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.x = element_text(face = "bold", size = 11),
    strip.text.y = element_text(face = "bold", size = 11, angle = 0),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "right",
    legend.text = element_text(size = 9), 
    legend.title = element_text(size = 10),
  )

print(p_stacked_depth)
ggsave("predicted_relative_abundance_stacked_bars_with_depth_facet.pdf",
       plot = p_stacked_depth, width = 14, height = 10)


### Fig. 4b - Absolute proportional error per tool ###

# reload necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)

#### reload data // reformat master2 if necessary ####
master <- read.csv("Simulated_metagenomes_all_references_raw_predicted_abundances.csv") # Available on FigShare

master2 <- master %>%
  mutate(
    sample_clean    = sub("^HiSeq_", "", sample_id),          # drop HiSeq_ prefix if present
    sample_with_rep = sample_clean,                           # keep replicate-aware sample id here
    sample_base     = sub("_rep[0-9]+$", "", sample_clean),  # drop _repN for plotting
    strain_category = case_when(
      strain_id == k12_id   ~ "K12-MG1655",
      strain_id == sakai_id ~ "O157:H7 Sakai",
      TRUE                  ~ "Other"
    )
  )

# 1. Restrict to Sakai AND K12 rows and valid detections
df_two <- master2 %>%
  filter(
    ecoli_strain %in% c("O157:H7 str.Sakai", "K-12 MG1655"),
    strain_present_in_sample == TRUE,
    detected == TRUE
  ) %>%
  mutate(
    absolute_proportional_error = as.numeric(absolute_proportional_error)
  ) %>%
  filter(!is.na(absolute_proportional_error))

# 2. First average triplicates per strain
df_sample_mean <- df_two %>%
  group_by(sample_base, strain_profiling_tool, ecoli_strain) %>%
  summarise(
    mean_APE_strain = mean(absolute_proportional_error),
    .groups = "drop"
  ) %>%
  group_by(sample_base, strain_profiling_tool) %>%
  summarise(
    mean_APE_sample = mean(mean_APE_strain),
    .groups = "drop"
  )

# 3. Compute mean ± SD per profiler
df_profiler_summary <- df_sample_mean %>%
  group_by(strain_profiling_tool) %>%
  summarise(
    mean_APE = mean(mean_APE_sample),
    sd_APE = sd(mean_APE_sample),
    n_samples = n(),
    .groups = "drop"
  )

print(df_profiler_summary)

# 4. Bar chart
p_mean_APE_bar <- ggplot(
  df_profiler_summary,
  aes(x = strain_profiling_tool, y = mean_APE)
) +
  geom_col(aes(fill = strain_profiling_tool), width = 0.7) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_APE - sd_APE),
      ymax = mean_APE + sd_APE),
    width = 0.2
  ) +
  geom_text(
    aes(
      y = mean_APE + sd_APE,
      label = sprintf("%.2f", mean_APE)),
    vjust = -0.6,
    size = 4
  ) +
  scale_y_continuous(
    name = "Mean absolute proportional error"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Strain-level profiling tool") +
  labs(
    x = "Strain-level profiling tool",
    #title = "Mean absolute proportional error (mean of K12 & Sakai)",
    #subtitle = "Mean ± SD across samples (triplicates averaged)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "",
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 37.5, hjust = 1),
    axis.title.x = element_text(vjust = -0.75)
  )

print(p_mean_APE_bar)

ggsave("mean_APE_bars_per_tool.pdf",
       plot = p_mean_APE_bar, width = 12, height = 6)


### cowplot for combined single pdf for Fig. 4 ###

# Fig.4
p4a <- p_stacked_depth
p4b <- p_mean_APE_bar

plot_grid(p4a, p4b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(2.5, 1))
ggsave("Simulated_metagenomes_Fig4.pdf", width = 16, height = 10)


#############################################################################

# Also part of this results narrative is Fig. S3 #
### Fig. S3a - Mapped reads PathoScope and Strainify ###

library(tidyverse)
library(readxl)

# Load files
all_refs_strainify_mapped_reads  <- read.csv("all_refs_strainify_reads_mapped.csv") # Both files available on FigShare.
all_refs_pathoscope_mapped_reads <- read.csv("all_refs_pathoscope_reads_mapped.csv")

# Combine into one dataframe
all_refs_mapped_reads_combined <- bind_rows(
  all_refs_strainify_mapped_reads,
  all_refs_pathoscope_mapped_reads
)

# Add replicate-stripped ID and depth
all_refs_mapped_reads_combined <- all_refs_mapped_reads_combined %>%
  mutate(
    replicate = str_extract(sample_id, "(?<=_rep)\\d+") %>% as.integer(),
    base_sample_id = str_remove(sample_id, "_rep\\d+$"),
    depth = str_extract(base_sample_id, "\\d+M")
  )

# Helper to convert sample_id into a compact ratio label
build_proportion_label <- function(sid) {
  sid_lower <- tolower(sid)
  
  # catch control / typo variants
  if (str_detect(sid_lower, "control|contorl")) {
    return("0.00:0.00")
  }
  
  sid <- str_remove(sid, "^HiSeq_")
  sid <- str_remove(sid, "_rep\\d+$")
  sid <- str_remove(sid, "^(20M|50M|100M)_")
  sid_clean <- str_remove_all(sid, "K12|Sakai|k12|sakai|K12_|Sakai_")
  
  nums <- str_extract_all(sid_clean, "\\d*\\.?\\d+")[[1]]
  if (length(nums) < 2) return(sid)
  
  p1 <- as.numeric(nums[length(nums) - 1])
  p2 <- as.numeric(nums[length(nums)])
  
  # convert percentages like 99 -> 0.99, but leave 0/1 alone
  if (!is.na(p1) && p1 > 1) p1 <- p1 / 100
  if (!is.na(p2) && p2 > 1) p2 <- p2 / 100
  
  # if both are present but do not sum to 1, infer the second value
  eps <- 1e-6
  if (!is.na(p1) && !is.na(p2) && abs((p1 + p2) - 1) > eps) {
    p2 <- 1 - p1
  }
  
  sprintf("%.2f:%.2f", p1, p2)
}

# Summarise mean percent mapped reads across technical replicates
all_refs_percent_mapped_reads_combined_summary <- all_refs_mapped_reads_combined %>%
  group_by(base_sample_id, strain_profiling_tool) %>%
  summarise(
    mean_percent_mapped = mean(percent_mapped_reads, na.rm = TRUE),
    sd_percent_mapped   = sd(percent_mapped_reads, na.rm = TRUE),
    n_rep               = sum(!is.na(percent_mapped_reads)),
    se_percent_mapped   = ifelse(n_rep > 0, sd_percent_mapped / sqrt(n_rep), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    ratio_label = vapply(base_sample_id, build_proportion_label, character(1)),
    ratio_value = as.numeric(str_extract(ratio_label, "^[0-9.]+"))
  ) %>%
  arrange(ratio_value) %>%
  mutate(
    ratio_label = factor(ratio_label, levels = unique(ratio_label))
  )

# Plot

p_all_refs_percent_mapped_reads <- ggplot(
  all_refs_percent_mapped_reads_combined_summary,
  aes(x = ratio_label, y = mean_percent_mapped, fill = strain_profiling_tool)
) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = sprintf("%.2f", mean_percent_mapped)),
    position = position_dodge(width = 0.8),
    vjust = -0.5, angle = 0, size = 3.5
  ) +
  geom_errorbar(
    aes(
      ymin = mean_percent_mapped - se_percent_mapped,
      ymax = mean_percent_mapped + se_percent_mapped
    ),
    position = position_dodge(width = 0.8),
    width = 0.3,
    na.rm = TRUE
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.title.x = element_text(vjust = 0.25),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9)
  ) +
  labs(
    x = "True K12-MG1655 : O157:H7 str. Sakai ratio in simulated metagenome sample",
    y = "Percentage mapped reads (within simulated metagenome)",
    fill = "Strain-level profiling tool"
  )

print(p_all_refs_percent_mapped_reads)
ggsave("mapped_reads_Strainify_Pathoscope.pdf",
       plot = p_all_refs_percent_mapped_reads, width = 12, height = 6)


### Fig. S3b - Mean absolute proportional error, faceted by true strains (K12/Sakai) ###

# reload necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(readxl)

#### reload data // reformat master2 if necessary ####
master <- read.csv("Simulated_metagenomes_all_references_raw_predicted_abundances.csv") # Available on FigShare

master2 <- master %>%
  mutate(
    sample_clean    = sub("^HiSeq_", "", sample_id),          # drop HiSeq_ prefix if present
    sample_with_rep = sample_clean,                           # keep replicate-aware sample id here
    sample_base     = sub("_rep[0-9]+$", "", sample_clean),  # drop _repN for plotting
    strain_category = case_when(
      strain_id == k12_id   ~ "K12-MG1655",
      strain_id == sakai_id ~ "O157:H7 Sakai",
      TRUE                  ~ "Other"
    )
  )

# 1. Restrict to Sakai AND K12 rows and valid detections
df_two <- master2 %>%
  filter(
    ecoli_strain %in% c("O157:H7 str.Sakai", "K-12 MG1655"),
    strain_present_in_sample == TRUE,
    detected == TRUE
  ) %>%
  mutate(
    absolute_proportional_error = as.numeric(absolute_proportional_error)
  ) %>%
  filter(!is.na(absolute_proportional_error))

# 2. First average triplicates per strain (KEEP strain separate)
df_sample_mean <- df_two %>%
  group_by(sample_base, strain_profiling_tool, ecoli_strain) %>%
  summarise(
    mean_APE_sample = mean(absolute_proportional_error),
    .groups = "drop"
  )

# 3. Compute mean ± SD per profiler *within each strain*
df_profiler_summary <- df_sample_mean %>%
  group_by(strain_profiling_tool, ecoli_strain) %>%
  summarise(
    mean_APE = mean(mean_APE_sample),
    sd_APE = sd(mean_APE_sample),
    n_samples = n(),
    .groups = "drop"
  )

print(df_profiler_summary)

# 4. Faceted bar chart (free y + labels above error bars)
p_mean_APE_bar_facet <- ggplot(
  df_profiler_summary,
  aes(x = strain_profiling_tool, y = mean_APE, fill = strain_profiling_tool)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_APE - sd_APE),
      ymax = mean_APE + sd_APE
    ),
    width = 0.2
  ) +
  geom_text(
    aes(
      y = mean_APE + sd_APE,
      label = sprintf("%.2f", mean_APE)
    ),
    vjust = -0.4,
    size = 4
  ) +
  facet_wrap(~ ecoli_strain, scales = "free_y") +
  scale_y_continuous(
    name = "Mean absolute proportional error",
    expand = expansion(mult = c(0, 0.1))  # adds headroom for labels
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Strain-level profiling tool") +
  labs(
    x = "Strain-level profiling tool"
    #title = "Mean absolute proportional error by strain",
    #subtitle = "Mean ± SD across samples (triplicates averaged); note different y-axis per strain"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(vjust = 0.25),
    axis.title.y = element_text(hjust = 0.25),
    legend.position = ""
  )

print(p_mean_APE_bar_facet)
ggsave("mean_APE_bars_per_tool_faceted_by_strain.pdf",
       plot = p_mean_APE_bar_facet, width = 12, height = 6)


### cowplot for combined single pdf for Fig. S3 ###

# Fig. S3
pS3a <- p_all_refs_percent_mapped_reads
pS3b <- p_mean_APE_bar_facet

plot_grid(pS3a, pS3b, labels = c("(a)", "(b)"), ncol = 2)
ggsave("Simulated_metagenomes_FigS3.pdf", width = 16, height = 8)