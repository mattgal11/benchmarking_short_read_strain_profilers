# ============================================================
# Fig. 9 Computational benchmarking of tools
#
# Metagenome benchmarking:
#   - 42 metagenomes per tool
#   - Boxplots show distribution across samples
#   - Individual points show each metagenome
#
# Database benchmarking:
#   - Single measurement per database-building step
#   - Bar plots used for runtime and peak memory

# Additional supplementary figures for this results section
#   - StrainR2 memory flag (-m) comparisons
#   - Paired metagenome file size distribution
#   - StrainScan output by metagenome size
# ============================================================

# ------------------------------------------------------------
# 1. Load necessary libraries
# ------------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(colorspace)
library(ggbreak)

# ============================================================
# METAGENOME COMPUTATIONAL BENCHMARKING
# ============================================================

# ------------------------------------------------------------
# 2. Load metagenome computational benchmarking data
# ------------------------------------------------------------
combined_comp_benchmark <- read.csv("real_world_spike_ins_comp_benchmark_profiling.csv") # <-- Available on FigShare.

# Sanity checks
head(combined_comp_benchmark)
table(combined_comp_benchmark$tool)

# ------------------------------------------------------------
# 3. Calculate averages of all numeric benchmark statistics
#    per tool
# ------------------------------------------------------------
comp_average_by_tool <- combined_comp_benchmark %>%
  mutate(
    user_s = as.numeric(user_s),
    system_s = as.numeric(system_s),
    cpu_pct = as.numeric(cpu_pct)
  ) %>%
  group_by(tool) %>%
  summarise(
    across(
      c(user_s,
        system_s,
        cpu_pct,
        wall_min,
        wall_h,
        peak_rss_kb,
        peak_rss_gib),
      ~ mean(.x, na.rm = TRUE)
    )
  )

print(comp_average_by_tool)

# Save table if required
# write.csv(comp_average_by_tool, "comp_averages_by_tool.csv", row.names = FALSE)

# ------------------------------------------------------------
# 4. Calculate summary statistics for key metrics
# Runtime = wall-clock time (hours)
# Memory  = peak RSS (GiB)
# ------------------------------------------------------------
key_comp_statistics_per_tool <- combined_comp_benchmark %>%
  group_by(tool) %>%
  summarise(
    n_wall = sum(!is.na(wall_h)),
    n_memory = sum(!is.na(peak_rss_gib)),
    # Wall-clock runtime
    mean_wall_h = mean(wall_h, na.rm = TRUE),
    sd_wall_h = sd(wall_h, na.rm = TRUE),
    median_wall_h = median(wall_h, na.rm = TRUE),
    Q1_wall_h = quantile(wall_h, 0.25, na.rm = TRUE),
    Q3_wall_h = quantile(wall_h, 0.75, na.rm = TRUE),
    IQR_wall_h = IQR(wall_h, na.rm = TRUE),
    min_wall_h = min(wall_h, na.rm = TRUE),
    max_wall_h = max(wall_h, na.rm = TRUE),
    # Peak memory
    mean_peak_rss_gib = mean(peak_rss_gib, na.rm = TRUE),
    sd_peak_rss_gib = sd(peak_rss_gib, na.rm = TRUE),
    median_peak_rss_gib = median(peak_rss_gib, na.rm = TRUE),
    Q1_peak_rss_gib = quantile(peak_rss_gib, 0.25, na.rm = TRUE),
    Q3_peak_rss_gib = quantile(peak_rss_gib, 0.75, na.rm = TRUE),
    IQR_peak_rss_gib = IQR(peak_rss_gib, na.rm = TRUE),
    min_peak_rss_gib = min(peak_rss_gib, na.rm = TRUE),
    max_peak_rss_gib = max(peak_rss_gib, na.rm = TRUE))

print(key_comp_statistics_per_tool)

# Save summary table if required
write.csv(key_comp_statistics_per_tool, "key_comp_statistics_per_tool_final.csv", row.names = FALSE)

# ------------------------------------------------------------
# 5. Plot wall-clock runtime across the 42 metagenomes
# Define colour for tools first (PanTax missing)
# Boxplot = median + IQR + whiskers
# Points  = individual metagenomes
# ------------------------------------------------------------
tool_colours <- c(
  "PanTax"                       = "#1B9E77",
  "PathoScope"                   = "#D95F02",
  "StrainGE"                     = "#7570B3",
  "Strainify"                    = "#E7298A",
  "StrainR2"                     = "#66A61E",
  "StrainScan"                   = "#E6AB02",
  "StrainScan (super low depth)" = "#A6761D")

comp_runtime_plot <- ggplot(combined_comp_benchmark,
    aes(x = tool, y = wall_h, fill = tool)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    linewidth = 0.5) +
  geom_point(
    aes(color = tool),
    position = position_jitter(
      width = 0.15, height = 0, seed = 123),
    size = 2, alpha = 0.7) +
  # Add median value as text
  stat_summary(
    fun = median, geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -1.75, fontface = "bold", size = 3.5, color = "black") +
  # Add mean as a white diamond (optional secondary statistic)
  #stat_summary(
  #  fun = mean, geom = "point",
  #  shape = 23, size = 3.5,
  #  fill = "white", color = "black") +
  scale_fill_manual(values = tool_colours) +
  scale_color_manual(values = tool_colours) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 15)}) +
  labs(
    x = "Strain-level profiling tool",
    y = "Wall-clock time (hours)",
    # title = expression("Profiling metagenomes (" * italic(n) * " = 42) wall clock time")) +
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9, face = "bold"))

print(comp_runtime_plot)

# Save if required
#ggsave("comp_runtime_boxplot.pdf", comp_runtime_plot, width = 8, height = 6)

# ------------------------------------------------------------
# 6. Plot peak memory across the 42 metagenomes
#
# Boxplot = median + IQR + whiskers
# Points  = individual metagenomes
# Diamond = arithmetic mean (optional)
# ------------------------------------------------------------
comp_memory_plot <- ggplot(combined_comp_benchmark,
  aes(x = tool, y = peak_rss_gib, fill = tool)) +
  geom_boxplot(
    width = 0.6, outlier.shape = NA, linewidth = 0.5) +
  geom_point(
    aes(color = tool),
    position = position_jitter(
      width = 0.15, height = 0, seed = 123),
    size = 2, alpha = 0.7) +
  # Add median value as text
  stat_summary(
    fun = median, geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.8, fontface = "bold", size = 3.5, color = "black") +
  # Add mean as a white diamond
  # stat_summary(
  #  fun = mean, geom = "point",
  #  shape = 23, size = 3.5,
  #  fill = "white", color = "black") +
  scale_fill_manual(values = tool_colours) +
  scale_color_manual(values = tool_colours) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 15)}) +
  scale_y_continuous(
    limits = c(0, 400),
    breaks = c(0, 5, 10, 15, 20, 25, 300, 350, 400),
    expand = expansion(mult = c(0, 0))) +
  ggbreak::scale_y_break(
    c(25, 300),
    scales = 0.15, space = 0.10, symbol = "slash") +
  labs(
    x = "Strain-level profiling tool",
    y = "Peak memory usage (GiB)",
    # title = expression("Profiling metagenomes (" * italic(n) * " = 42) peak memory usage")) +
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9, face = "bold"))

print(comp_memory_plot)

# Save if required
#ggsave("comp_memory_boxplot.pdf", comp_memory_plot, width = 8, height = 6)


# ============================================================
# DATABASE-BUILD COMPUTATIONAL BENCHMARKING
# ============================================================

# ------------------------------------------------------------
# 7. Load DB-building computational benchmarking data
# ------------------------------------------------------------
DB_comp_benchmark <- read.csv("real_world_spike_ins_comp_benchmark_DB_step.csv") # <-- Available on FigShare.

head(DB_comp_benchmark)

# ------------------------------------------------------------
# 8. Create separate plotting category for PanTax index step
# Also have a combined label
# ------------------------------------------------------------
DB_comp_benchmark <- DB_comp_benchmark %>%
  mutate(
    plot_group = ifelse(
      tool == "PanTax" & step == "DB_index",
      "PanTax index", tool),
    # Control stacking order for runtime plot
    plot_group = factor(
      plot_group,
      levels = c(
        "PanTax index",
        "PanTax",
        "PathoScope",
        "StrainGE",
        "Strainify",
        "StrainR2",
        "StrainScan",
        "StrainScan (memory efficient)")))

pantax_total <- DB_comp_benchmark %>%
  filter(tool == "PanTax") %>%
  summarise(
    tool = "PanTax",
    wall_h = sum(wall_h, na.rm = TRUE))

# ------------------------------------------------------------
# 9. (Re-)Define tool colours
# One Dark2 colour per tool.
# PanTax index step uses lighter PanTax colour.
# ------------------------------------------------------------
dark2 <- brewer.pal(7, "Dark2")

tool_colours <- c(
  "PanTax"= dark2[1],
  "PanTax index"= lighten(dark2[1], 0.4),
  "PathoScope"= dark2[2],
  "StrainGE"= dark2[3],
  "Strainify" = dark2[4],
  "StrainR2"= dark2[5],
  "StrainScan"= dark2[6],
  "StrainScan (memory efficient)" = dark2[7])

# ------------------------------------------------------------
# 10. Plot DB-building wall-clock runtime
# PanTax DB build + index are stacked because runtime is
# additive across the sequential DB-building stages.
# ------------------------------------------------------------
DB_comp_runtime_plot <- ggplot(DB_comp_benchmark,
  aes(x = tool, y = wall_h, fill = plot_group)) +
  geom_col() +
   geom_text(
     data = DB_comp_benchmark %>%
       filter(tool == "PanTax"),
     aes(label = round(wall_h, 2)),
     position = position_stack(vjust = 0.5),
     size = 3.5, colour = "white", fontface = "bold") +
  geom_text(data = pantax_total,
    aes(x = tool, y = wall_h,
        label = round(wall_h, 2)),
    vjust = -0.5, size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = DB_comp_benchmark %>%
      filter(tool != "PanTax"),
    aes(label = round(wall_h, 2),
      y = wall_h),
    vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = tool_colours) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 15)}) +
  labs(
    x = "",
    y = "Wall-clock time (hours)",
    #title = "Database-building step wall-clock time") +
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 32.5, hjust = 0.5, size = 9, face = "bold"))

print(DB_comp_runtime_plot)

# Save if required
# ggsave("DB_step_comp_runtime_plot.png", DB_comp_runtime_plot, width = 6, height = 4, dpi = 300)

# ------------------------------------------------------------
# 11. Define separate PanTax order for memory plot
# Peak RSS is NOT additive, so PanTax DB build and index
# are plotted as neighbouring bars rather than stacked.
# ------------------------------------------------------------
DB_comp_benchmark$plot_group_memory <- factor(
  DB_comp_benchmark$plot_group,
  levels = c(
    "PanTax",
    "PanTax index",
    "PathoScope",
    "StrainGE",
    "Strainify",
    "StrainR2",
    "StrainScan",
    "StrainScan (memory efficient)")
  )

# ------------------------------------------------------------
# 12. Plot DB-building peak memory
# ------------------------------------------------------------
DB_comp_memory_plot <- ggplot(DB_comp_benchmark,
  aes(x = tool, y = peak_rss_gib,
    fill = plot_group_memory)) +
  geom_col(position = position_dodge(width = 0.8)) +
  # Labels above each bar
  geom_text(
    aes(label = round(peak_rss_gib, 2),
      y = peak_rss_gib), position = position_dodge(width = 0.8),
    vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = tool_colours) +
  scale_x_discrete(
    labels = function(x) {
      stringr::str_wrap(x, width = 15)}) +
  labs(
    x = "",
    y = "Peak memory usage (GiB)",
    # title = "Database-building step peak memory usage") +
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x = element_text(angle = 32.5, hjust = 1, size = 9, face = "bold"))

print(DB_comp_memory_plot)

# Save if required
# ggsave("DB_step_comp_memory_plot.png", DB_comp_memory_plot, width = 6, height = 4, dpi = 300)


### COWPLOT FOR COMBINED FIG8 (4 PANELS) ###

library(cowplot)
library(grid)

p9a <- DB_comp_runtime_plot
p9b <- DB_comp_memory_plot
p9c <- comp_runtime_plot
p9d <- grid.grabExpr(
  print(comp_memory_plot)) # Needs special formatting to preserve the y axis break

combined_fig9 <- plot_grid(
  p9a, p9b, p9c, p9d,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  ncol = 2, rel_widths = c(1, 1),
  label_x = -0.01, label_y = 1.01)

print(combined_fig9)
ggsave("combined_Fig9.pdf", combined_fig9, width = 16, height = 10)


# ============================================================
# ADDITIONAL STRAINR2 PLOT ACROSS TESTED MEMORY (-m) FLAGS
# ============================================================

# StrainR2 specifically allows control of -m flag (BBMap step) which influences peak RSS being measured here.
# We ran strainR2 across different -m flags to see how this influences peak memory usage across the run.
# This is an additional supplementary figure.
# Note strainR2 failed with -m 8.

# ------------------------------------------------------------
# 1. Load data and calculate averages per -m flag
# ------------------------------------------------------------
strainR2_memory_comp <- read.csv("strainR2_memory_comparisons.csv") # <-- Available on FigShare.

head(strainR2_memory_comp)
table(strainR2_memory_comp$tool)

# Calculate averages per strainR2 memory flag
strainR2_average_per_mem_flag <- strainR2_memory_comp %>%
  group_by(tool) %>%
  summarise(
    across(
      c(user_s, system_s, cpu_pct, wall_min, wall_h, peak_rss_kb, peak_rss_gib),
      ~ mean(.x, na.rm = TRUE)))

# Save table of strainR2 -m flag averages per tool to CSV
# write.csv(strainR2_average_per_mem_flag, "strainR2_average_per_mem_flag.csv", row.names = FALSE)

# ------------------------------------------------------------
# 2. Calculate mean and s.d. for wall_h and peak_rss_gib
# ------------------------------------------------------------
strainR2_key_comp_statistics_per_mem_flag <- strainR2_memory_comp %>%
  group_by(tool) %>%
  summarise(
    n = n(),
    mean_wall_h = mean(wall_h, na.rm = TRUE),
    sd_wall_h = sd(wall_h, na.rm = TRUE),
    mean_peak_rss_gib = mean(peak_rss_gib, na.rm = TRUE),
    sd_peak_rss_gib = sd(peak_rss_gib, na.rm = TRUE))

# ------------------------------------------------------------
# 3. Plot wall-clock runtime across StrainR2 -m flags
# As box plots to show the distribution across 42 samples
# ------------------------------------------------------------
comp_runtime_strainR2_mem_flags <- ggplot(strainR2_memory_comp,
  aes(x = tool, y = wall_h, fill = tool)) +
  geom_boxplot(width = 0.6, outlier.shape = NA,linewidth = 0.5) +
  geom_point(aes(color = tool),
    position = position_jitter(
      width = 0.15, height = 0, seed = 123),
    size = 2, alpha = 0.7) +
  # Median value
  stat_summary(
    fun = median, geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.8, fontface = "bold",
    size = 3.5, color = "black") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "StrainR2 -m flag (BBMap step)",
       y = "Wall-clock time (hours)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x = element_text(size = 10, face = "bold"))

print(comp_runtime_strainR2_mem_flags)
ggsave("strainR2_runtime_by_mem_flag_box_plot.pdf", comp_runtime_strainR2_mem_flags, width = 8, height = 6)

# ------------------------------------------------------------
# 4. Plot peak memory across StrainR2 -m flags
# ------------------------------------------------------------
comp_memory_strainR2_mem_flags <- ggplot(strainR2_memory_comp,
  aes(x = tool, y = peak_rss_gib, fill = tool)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.5) +
  geom_point(aes(color = tool),
    position = position_jitter(
      width = 0.15, height = 0, seed = 123),
    size = 2, alpha = 0.7) +
  # Median value
  stat_summary(
    fun = median, geom = "text",
    aes(label = after_stat(sprintf("%.2f", y))),
    vjust = -0.8, fontface = "bold",
    size = 3.5, color = "black") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "StrainR2 -m flag (BBMap step)",
       y = "Peak memory usage (GiB)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.text.x = element_text(size = 10, face = "bold"))

print(comp_memory_strainR2_mem_flags)
ggsave("strainR2_memory_by_mem_flag_box_plot.pdf", comp_memory_strainR2_mem_flags, width = 8, height = 6)

# ------------------------------------------------------------
# 5. Save as cowplot for supplementary figure
# ------------------------------------------------------------
library(cowplot)

S7a <- comp_runtime_strainR2_mem_flags
S7b <- comp_memory_strainR2_mem_flags

plot_grid(S7a, S7b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(1, 1))
ggsave("combined_FigS7.pdf", width = 12, height = 8)


# ============================================================
# ADDITIONAL plot: Paired metagenome file size (R1 + R2)
# ============================================================

# ------------------------------------------------------------
# 1. (re)Load libraries and data
# ------------------------------------------------------------
library(readxl)
library(tidyverse)

metagenome_file_sizes <- read.csv("real_world_spike_ins_metagenome_file_sizes.csv") # <-- Available on FigShare.

# ------------------------------------------------------------
# 2. Calculate total paired metagenome size (R1 + R2) and summary statistics
# ------------------------------------------------------------
paired_metagenome_sizes <- metagenome_file_sizes %>%
  mutate(
    sample = str_remove(Filename, "_R[12]\\.fastq\\.gz$")
  ) %>%
  group_by(sample) %>%
  summarise(
    paired_size_gb = sum(`Size (GB)`),
    .groups = "drop")

# Summary statistics
mean_size <- mean(paired_metagenome_sizes$paired_size_gb)
min_size  <- min(paired_metagenome_sizes$paired_size_gb)
max_size  <- max(paired_metagenome_sizes$paired_size_gb)

cat("Mean paired metagenome size:", round(mean_size, 2), "Gb\n")
cat("Range:",
    round(min_size, 2), "-",
    round(max_size, 2), "Gb\n")

# ------------------------------------------------------------
# 3. Plot paired metagenome sizes as whisker box plot
# ------------------------------------------------------------
p_paired_metagenome_sizes <- ggplot(paired_metagenome_sizes, aes(x = "", y = paired_size_gb)) +
  geom_boxplot(width = 0.25, outlier.shape = NA,
               fill = "#BFD7EA", colour = "#2F4858", linewidth = 0.8) +
  geom_jitter(width = 0.06, size = 2, alpha = 0.65, colour = "#2F4858") +
  geom_hline(yintercept = mean_size, linetype = "dashed") +
  annotate("text", x = 1.22, y = 2.8,
           label = paste0("Mean = ", round(mean_size, 2), " GB"), hjust = 0, fontface = "bold") +
  labs(x = NULL, y = "Paired metagenome size (GB)",
       title = "Paired metagenome file size",
       subtitle = paste0(
         "R1 + R2 after quality control | Range = ",
         round(min_size, 2), "-",
         round(max_size, 2), " GB")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Print and save plot
print(p_paired_metagenome_sizes)
ggsave("FigS6.pdf", p_paired_metagenome_sizes, width = 8, height = 6)


# ============================================================
# StrainScan sub-analysis:
# Are samples where only one strain was detected associated
# with smaller metagenome size / lower sequencing depth?
# ============================================================

# ------------------------------------------------------------
# 1. Load libraries and data
# ------------------------------------------------------------
library(readxl)
library(tidyverse)

master <- read.csv("real_world_spike_ins_benchmarking_results.csv") # <-- Available on FigShare.
metagenome_file_sizes <- read.csv("real_world_spike_ins_metagenome_file_sizes.csv") # <-- Available on FigShare.

# ------------------------------------------------------------
# 2. (re)Calculate total paired metagenome size (R1 + R2)
# ------------------------------------------------------------
paired_sizes <- metagenome_file_sizes %>%
  mutate(
    sample_id = str_remove(
      Filename,
      "_final_spiked_R[12]\\.fastq\\.gz$"
    )
  ) %>%
  group_by(sample_id) %>%
  summarise(
    paired_size_gb = sum(`Size (GB)`),
    .groups = "drop")

# ------------------------------------------------------------
# 3. Count how many strains StrainScan detected per sample
# ------------------------------------------------------------
strainscan <- master %>%
  filter(
    strain_profiling_tool == "StrainScan",
    detected == TRUE
  ) %>%
  count(
    sample_id,
    name = "n_strains"
  ) %>%
  mutate(
    result = if_else(
      n_strains == 1,
      "One strain detected",
      "Multiple strains detected"))

# ------------------------------------------------------------
# 4. Join StrainScan results with metagenome size
# ------------------------------------------------------------
strainscan_size <- strainscan %>%
  left_join(
    paired_sizes,
    by = "sample_id")

# Set plotting order
strainscan_size <- strainscan_size %>%
  mutate(
    result = factor(
      result,
      levels = c(
        "Multiple strains detected",
        "One strain detected")))

# ------------------------------------------------------------
# 5. Summary statistics
# ------------------------------------------------------------
summary_stats <- strainscan_size %>%
  group_by(result) %>%
  summarise(
    n = n(),
    mean_GB = mean(paired_size_gb, na.rm = TRUE),
    median_GB = median(paired_size_gb, na.rm = TRUE),
    min_GB = min(paired_size_gb, na.rm = TRUE),
    max_GB = max(paired_size_gb, na.rm = TRUE),
    .groups = "drop")

print(summary_stats)

# ------------------------------------------------------------
# 6. Wilcoxon rank-sum test
# ------------------------------------------------------------
wilcox_result <- wilcox.test(
  paired_size_gb ~ result,
  data = strainscan_size,
  exact = FALSE)

print(wilcox_result)
p_value <- wilcox_result$p.value

# ------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------
p_strainscan_size <- ggplot(
  strainscan_size,
  aes(x = result, y = paired_size_gb)) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "#BFD7EA",
    colour = "#2F4858", linewidth = 0.8) +
  geom_jitter(width = 0.06, size = 2, alpha = 0.65, colour = "#2F4858") +
  annotate("text",
    x = 1.5,
    y = max(strainscan_size$paired_size_gb, na.rm = TRUE) * 1.06,
    label = paste0(
      "Wilcoxon p = ",
      format.pval(p_value, digits = 3)),
    fontface = "bold") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL,
       y = "Paired metagenome size (GB)",
       title = "StrainScan output by metagenome size",
       subtitle = "R1 + R2 metagenome file size") +
  theme_minimal(base_size = 12) +
  theme(axis.title.y = element_text(
      margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank())

# ------------------------------------------------------------
# 8. Print and save plot
# ------------------------------------------------------------
print(p_strainscan_size)
ggsave("Fig.S4.pdf", p_strainscan_size, width = 8, height = 6)