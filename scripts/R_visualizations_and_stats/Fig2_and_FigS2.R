# Final plotting scripts for the 99% human 1% ecoli dataset plots included in the draft manuscript (Fig. 2 and Fig. S2)

### Fig. 2a - Raw predicted abundance values per base tool (collapsed) ###

library(readxl)
library(dplyr)
library(ggplot2)

# Load the dataset
human99_ecoli1 <- read.csv("SRR13355226_raw_predicted_abundance_with_and_without_human_read_removal.csv") # Dataset available on FigShare. 

# -------- Rename and coerce types --------
human99_ecoli1 <- human99_ecoli1 %>%
  rename(
    tool              = Strain_profiling_tool_and_version,
    base_tool         = Base_tool,
    strain            = Escherichia_coli_strain,
    depth_by_coverage = Strain_depth_by_coverage.x.,
    error             = error,
    abs_err           = Absolute_error,
    truth             = Relative_ground_truth_abundance,
    ape               = Absolute_proportional_error
  ) %>%
  mutate(
    error = as.numeric(error),
    abs_err = as.numeric(abs_err),
    truth = as.numeric(truth),
    ape = as.numeric(ape),
    Predicted_abundance = as.numeric(Predicted_abundance),
    depth_by_coverage = as.numeric(as.character(depth_by_coverage)),
    strain_label = paste0(strain, " (", depth_by_coverage, "x)")
  ) %>%
  mutate(
    strain_label = factor(
      strain_label,
      levels = unique(strain_label[order(-depth_by_coverage)])
    )
  )

# Remove unwanted rows
human99_ecoli1v2 <- human99_ecoli1[-c(33:36, 41:44), ]

# -------- Collapse across conditions --------
collapsed_df <- human99_ecoli1v2 %>%
  group_by(base_tool, strain_label, depth_by_coverage) %>%
  summarise(
    Predicted_abundance = mean(Predicted_abundance, na.rm = TRUE),
    detected = ifelse(any(detected == "No"), "No", "Yes"),
    .groups = "drop"
  )

# -------- Plot (x = base_tool) --------

dodge <- position_dodge(width = 0.9)

max_y <- max(collapsed_df$Predicted_abundance, na.rm = TRUE)
if (!is.finite(max_y) || max_y <= 0) {
  upper_limit <- 1
  offset <- 0.05
} else {
  upper_limit <- max_y * 1.15
  offset <- max_y * 0.05
}

p_human99_pred_abundance_collapsed <- ggplot(
  collapsed_df,
  aes(x = base_tool, y = Predicted_abundance, fill = strain_label)
  ) +
  geom_col(position = dodge, width = 0.8, colour = NA) +
  
  # numeric labels
  geom_text(
    aes(label = sprintf("%.3f", Predicted_abundance)),
    position = dodge,
    hjust = -0.2,
    fontface = "bold",
    size = 3.0,
    angle = 90,
    na.rm = TRUE
  ) +
  
  # red stars
  geom_text(
    aes(
      label = ifelse(detected == "No", "*", NA),
      y = ifelse(detected == "No",
                 Predicted_abundance + offset,
                 NA)
    ),
    position = dodge,
    colour = "#B22222",
    size = 6,
    vjust = -0.25,
    na.rm = TRUE
  ) +
  
  scale_fill_viridis_d(
    name = expression(italic(E.~coli)~" strain (sequencing depth (x), ground truth)"),
    option = "turbo"
  ) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(
    x = "Strain-level profiling tool",
    y = "Predicted relative abundance"
  ) +
  scale_y_continuous(
    limits = c(0, upper_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text  = element_text(size = 9),
    axis.text.x = element_text(angle = 22.5, size = 9, vjust = 0.5, hjust = 0.5),
    plot.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(t = 10, r = 20, b = 30, l = 10)
  )

print(p_human99_pred_abundance_collapsed)

ggsave(
  "99human_1ecoli_figures/99human_1ecoli_pred_abundance_collapsed_base_tool.pdf",
  plot = p_human99_pred_abundance_collapsed, width = 14, height = 10)


### Fig. 2b - mean APE per tool ###

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)

# Reload and clean data
human99_ecoli1 <- read.csv("SRR13355226_raw_predicted_abundance_with_and_without_human_read_removal.csv") %>%
  rename(
    tool              = Strain_profiling_tool_and_version,
    base_tool         = Base_tool,
    strain            = Escherichia_coli_strain,
    depth_by_coverage = Strain_depth_by_coverage.x.,
    error             = error,
    abs_err           = Absolute_error,
    truth             = Relative_ground_truth_abundance,
    ape               = Absolute_proportional_error
  ) %>%
  mutate(
    error = as.numeric(error),
    abs_err = as.numeric(abs_err),
    truth = as.numeric(truth),
    ape = as.numeric(ape),
    depth_by_coverage = as.numeric(as.character(depth_by_coverage))
  ) %>%
  filter(
    !is.na(base_tool),
    base_tool != "Ground truth"
  )

# remove unwanted rows
human99_ecoli1 <- human99_ecoli1[-c(33:36, 41:44), ]

# Summarise mean ± SD
summary_base_ape <- human99_ecoli1 %>%
  group_by(base_tool) %>%
  summarise(
    n_used = sum(!is.na(ape) & ape != 0),
    n_not_detected = sum(detected == "No", na.rm = TRUE),
    mean_ape = ifelse(
      n_used > 0,
      mean(ape[!is.na(ape) & ape != 0]),
      NA_real_
    ),
    sd_ape = ifelse(
      n_used > 1,
      sd(ape[!is.na(ape) & ape != 0]),
      NA_real_
    ),
    sem_ape = ifelse(
      n_used > 1,
      sd_ape / sqrt(n_used),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    mean_plot = ifelse(is.na(mean_ape), 0, mean_ape),
    ymin = pmax(mean_plot - sd_ape, 0, na.rm = TRUE),
    ymax = mean_plot + sd_ape
  )

# Dynamic spacing
overall_max_mean <- max(summary_base_ape$mean_plot, na.rm = TRUE)
if (!is.finite(overall_max_mean) || overall_max_mean <= 0) overall_max_mean <- 1

label_offset <- overall_max_mean * 0.05
star_offset  <- overall_max_mean * 0.12

summary_base_ape <- summary_base_ape %>%
  mutate(
    mean_label = ifelse(
      !is.na(mean_ape),
      number(mean_ape, accuracy = 0.01),
      "NA"
    ),
    label_y = ymax + label_offset,
    star = ifelse(n_not_detected > 0, "*", NA),
    star_y = ymax + star_offset
  )

# Plot
p_mean_ape_base <- ggplot(
  summary_base_ape,
  aes(x = base_tool, y = mean_plot, fill = base_tool)
  ) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(
    aes(ymin = ymin, ymax = ymax),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = mean_label, y = label_y),
    vjust = 0,
    size = 3.75,
    fontface = "bold"
  ) +
  geom_text(
    aes(label = star, y = star_y),
    colour = "#B22222",
    size = 6,
    vjust = -0.5,
    na.rm = TRUE
  ) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "Strain-level profiling tool",
    y = "Mean Absolute Proportional Error (± SD)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 22.5, hjust = 0.75),
    legend.position = "none"
  )

print(p_mean_ape_base)

ggsave("99human_1ecoli_figures/mean_APE_base_tools.pdf",
       plot = p_mean_ape_base, width = 12, height = 6)


### cowplot joint multi-panel Fig.2 for final publication ###

library(cowplot)

# Fig.2
p2a <- p_human99_pred_abundance_collapsed
p2b <- p_mean_ape_base

plot_grid(p2a, p2b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(2, 1), label_x = c(0, -0.04))
ggsave("99human_1ecoli_figures/99human1ecoli_Fig2_combined.pdf", width = 16, height = 8)



##########################################################################################################

### Fig. S2a - APE per tool per strain (faceted); included here as it is part of this results narraitve ###

# reload necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)

# set tool order to match those above
base_tool_levels <- c(
  "PanTax",
  "PathoScope",
  "StrainGE",
  "Strainify",
  "StrainR2",
  "StrainScan",
  "StrainScan super low depth"
)

# reload the dataset and clean names
human99_ecoli1 <- read.csv("SRR13355226_raw_predicted_abundance_with_and_without_human_read_removal.csv") %>%
  rename(
    tool              = Strain_profiling_tool_and_version,
    base_tool         = Base_tool,
    strain            = Escherichia_coli_strain,
    depth_by_coverage = Strain_depth_by_coverage.x.,
    error             = error,
    abs_err           = Absolute_error,
    truth             = Relative_ground_truth_abundance,
    ape               = Absolute_proportional_error
  ) %>%
  mutate(
    error = as.numeric(error),
    abs_err = as.numeric(abs_err),
    truth = as.numeric(truth),
    ape = as.numeric(ape),
    depth_by_coverage = as.numeric(as.character(depth_by_coverage)),
    strain_label = paste0(strain, " (", depth_by_coverage, "x)"),
    base_tool = factor(base_tool, levels = base_tool_levels)
  ) %>%
  filter(!is.na(base_tool))

# remove data not of interest and compute missing APE if needed
human99_ecoli1 <- human99_ecoli1[-c(33:36, 41:44, 65:68), ]
human99_ecoli1 <- human99_ecoli1 %>%
  mutate(
    ape = ifelse(
      is.na(ape) & !is.na(abs_err) & !is.na(truth) & truth != 0,
      abs_err / truth,
      ape
    )
  )

# define strain label and order strains by coverage for plotting
strain_levels <- human99_ecoli1 %>%
  distinct(strain, depth_by_coverage, strain_label) %>%
  arrange(desc(depth_by_coverage)) %>%
  pull(strain_label)

plot_df <- human99_ecoli1 %>%
  mutate(
    strain_label = factor(strain_label, levels = strain_levels),
    base_tool = factor(base_tool, levels = base_tool_levels)
  ) %>%
  group_by(strain_label, strain, depth_by_coverage, base_tool) %>%
  summarise(
    ape = mean(ape, na.rm = TRUE),
    abs_err = mean(abs_err, na.rm = TRUE),
    truth = mean(truth, na.rm = TRUE),
    detected = ifelse(any(detected == "No", na.rm = TRUE), "No", "Yes"),
    .groups = "drop"
  ) %>%
  mutate(base_tool = factor(base_tool, levels = base_tool_levels))

per_facet <- plot_df %>%
  group_by(strain_label) %>%
  summarise(facet_max = max(ape, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    facet_max = ifelse(!is.finite(facet_max) | facet_max <= 0, 0.05, facet_max),
    facet_offset = pmax(facet_max * 0.05, 0.01),
    facet_ylim_top = facet_max + facet_offset * 3
  )

plot_df <- plot_df %>%
  left_join(per_facet, by = "strain_label") %>%
  mutate(
    star_y = ifelse(detected == "No", ape + facet_offset, NA_real_),
    label_y = ape
  )

stars_df <- plot_df %>%
  filter(!is.na(star_y) & is.finite(star_y))

# plotting
pd <- position_dodge2(width = 0.8, preserve = "single")
p_human99_ecoli1_base_tools_strain_facet <- ggplot(
  plot_df,
  aes(x = base_tool, y = ape, fill = base_tool, group = base_tool)
  ) +
  geom_col(width = 0.7, show.legend = FALSE, na.rm = TRUE, position = pd) +
  geom_text(
    aes(label = round(ape, 2), y = label_y),
    position = pd,
    vjust = -0.5,
    size = 3.75,
    fontface = "bold",
    na.rm = TRUE
  ) +
  geom_text(
    data = stars_df,
    aes(x = base_tool, y = star_y, label = "*", group = base_tool),
    position = pd,
    vjust = -0.5,
    colour = "#B22222",
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_blank(aes(y = facet_ylim_top)) +
  facet_wrap(~ strain_label, scales = "free_y", ncol = 2) +
  scale_x_discrete(limits = base_tool_levels, drop = FALSE) +
  labs(
    x = "Strain-level profiling tool",
    y = "Absolute Proportional Error"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none",
    panel.spacing = unit(0.6, "lines"),
    panel.background = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "white", colour = NA)
  )

# print and save
print(p_human99_ecoli1_base_tools_strain_facet)
ggsave("99human_1ecoli_figures/APE_per_base_tool_strain_facet.pdf",
  plot = p_human99_ecoli1_base_tools_strain_facet, width = 12, height = 6)


### Fig. S2b - Sensitivity per base tool ###

# Plots sensitivity of tools in 99% human 1% ecoli dataset

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Reload data and remove rows not of interest
human99_ecoli1 <- read.csv("SRR13355226_raw_predicted_abundance_with_and_without_human_read_removal.csv")
human99_ecoli1 <- human99_ecoli1[-c(33:36, 41:44, 65:68), ]

# count TP and FN per tool
sens_table_99human_1ecoli <- human99_ecoli1 %>%
  count(Strain_profiling_tool_and_version, call.type) %>%
  filter(call.type %in% c("TP", "FN")) %>%
  pivot_wider(names_from = call.type,
              values_from = n,
              values_fill = 0) %>%
  mutate(
    sensitivity = TP / (TP + FN))

# count TP and FN per tool
sens_table_base_tools_99human_1ecoli <- human99_ecoli1 %>%
  count(Base_tool, call.type) %>%
  filter(call.type %in% c("TP", "FN")) %>%
  pivot_wider(names_from = call.type,
              values_from = n,
              values_fill = 0) %>%
  mutate(
    sensitivity = TP / (TP + FN))

# plot sensitivity
p_sensitivity_base_tools_human99_1ecoli <- ggplot(sens_table_base_tools_99human_1ecoli,
                                                  aes(x = Base_tool, y = sensitivity, fill = Base_tool)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  #geom_text(aes(label = sensitivity),
  #          vjust = -0.5,
  #          size = 3.5,
  #          fontface = "bold") +
  scale_y_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(
    x = "Strain-level profiling tool",
    y = "Sensitivity"
    #title = "Sensitivity per tool in SPX9716028 (99% human, 1% E. coli)",
    #subtitle = "Sensitivity calculated as TP / (TP + FN)"
  ) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

print(p_sensitivity_base_tools_human99_1ecoli)
ggsave("99human_1ecoli_figures/sensitivty_per_base_tool.pdf",
       plot = p_sensitivity_base_tools_human99_1ecoli, width = 12, height = 6)


### Fig. S2c - Paired Wilcox stats test on absolute error ###

# load libraries
library(readxl)
library(dplyr)
library(tidyr)

# update df
stats_99human_1ecoli <- read.csv("SRR13355226_raw_predicted_abundance_with_and_without_human_read_removal.csv")
stats_99human_1ecoli <- stats_99human_1ecoli[-c(33:36, 41:44, 65:68), ]

df_wide <- stats_99human_1ecoli %>%
  select(Base_tool,
         Escherichia_coli_strain,
         Human_reads_removed,
         Absolute_error,
         Absolute_proportional_error,
         Predicted_abundance) %>%
  pivot_wider(
    names_from = Human_reads_removed,
    values_from = c(Absolute_error,
                    Absolute_proportional_error,
                    Predicted_abundance)
  )

# stats test: per tool
results <- df_wide %>%
  group_by(Base_tool) %>%
  summarise(
    mean_diff = mean(Absolute_error_Yes - Absolute_error_No),
    all_zero  = all(Absolute_error_Yes - Absolute_error_No == 0, na.rm = TRUE),
    p_value = if(all_zero) {
      NA_real_
    } else {
      wilcox.test(Absolute_error_Yes,
                  Absolute_error_No,
                  paired = TRUE)$p.value
    }
  )

results

# plot

library(dplyr)
library(ggplot2)

# Per-tool p-values (3 d.p.)
pvals <- df_wide %>%
  group_by(Base_tool) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(Absolute_error_Yes,
                  Absolute_error_No,
                  paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(
      is.na(p_value),
      "Wilcoxon p = NA\n(identical values)",
      paste0("Wilcoxon p = ",
             format(round(p_value, 3), nsmall = 3))
    )
  )

p_absolute_error <- ggplot(
  df_wide,
  aes(x = Absolute_error_No,
      y = Absolute_error_Yes,
      color = Base_tool)
) +
  geom_point(size = 3, alpha = 0.9, show.legend = FALSE) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed",
              colour = "black") +
  coord_equal(clip = "off") +
  facet_wrap(~ Base_tool, ncol = 7,
             labeller = labeller(Base_tool = label_wrap_gen(width = 20))) +
  geom_text(
    data = pvals,
    aes(label = label),
    x = Inf,
    y = -Inf,
    hjust = 1.05,
    vjust = -0.5,
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = "Absolute error (without human read removal)",
    y = "Absolute error (with human read removal)"
    #title = "Paired comparison of absolute error by tool",
    #subtitle = "Dashed line represents the identity line (y = x)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    plot.margin = margin(10, 20, 10, 10),
    axis.text.x = element_text(size = 7.75),
    axis.text.y = element_text(size = 7.75),
    axis.title.x = element_text(size = 12, vjust = -1),
    axis.title.y = element_text(size = 10)
  )

print(p_absolute_error)

ggsave("99human_1ecoli_figures/absolute_error_comparison_all_tools.pdf", 
       plot = p_absolute_error, width = 12, height = 6)

# check magnitude of differences
df_wide %>%
  mutate(diff = Absolute_error_Yes - Absolute_error_No) %>%
  summarise(
    median_diff = median(diff),
    max_diff = max(abs(diff))
  )

### cowplot joint multi-panel Fig. S2 for final publication ###

# Fig.S2
pS2a <- p_absolute_error
pS2b <- p_sensitivity_base_tools_human99_1ecoli
pS2c <- p_human99_ecoli1_base_tools_strain_facet

top_row <- plot_grid(pS2a, labels = c("(a)"), ncol = 1, label_y = 0.96)
bottom_row <- plot_grid(pS2b, pS2c, labels = c("(b)", "(c)"), ncol = 2, rel_widths = c(1, 2), label_y = c(1.06, 1.06))

# Combine and save
plot_grid(top_row, bottom_row, ncol = 1)
ggsave("99human_1ecoli_figures/99human1ecoli_FigS2_combined.pdf", width = 16, height = 8)