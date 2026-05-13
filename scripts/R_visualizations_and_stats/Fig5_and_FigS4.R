# Final plotting scripts for Fig. 5 and Fig. S4 (same results narrative)

### Fig. 5. Stacked bars by strain (phylogroup); Dataset: K12 & Sakai removed from reference database

library(dplyr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(scales)

# ---- Helper: build compact proportion label from sample id
build_proportion_label <- function(sid) {
  depth <- str_extract(sid, "20M|50M|100M")
  s_no_depth <- if (!is.na(depth)) str_remove(sid, paste0("^", depth, "_")) else sid
  
  s_clean <- str_remove_all(s_no_depth, "K12|Sakai|k12|sakai|K12_|Sakai_")
  nums <- str_extract_all(s_clean, "\\d*\\.?\\d+")[[1]]
  
  if (length(nums) < 2) return(s_no_depth)
  
  p1 <- as.numeric(nums[length(nums) - 1])
  p2 <- as.numeric(nums[length(nums)])
  
  if (!is.na(p1) && !is.na(p2) && (p1 > 1 || p2 > 1)) {
    total <- p1 + p2
    if (total > 0) {
      p1 <- p1 / total
      p2 <- p2 / total
    }
  }
  
  sprintf("%.2f:%.2f", p1, p2)
}

# Load and prepare data
master_no_K12_sakai <- read.csv("Simulated_metagenomes_K12_and_Sakai_removed_from_reference_database_raw_predicted_abundance.csv") # Available on FigShare.

master2_no_K12_sakai <- master_no_K12_sakai %>%
  mutate(
    sample_clean = sub("^HiSeq_", "", sample_id),
    sample_base  = sub("_rep[0-9]+$", "", sample_clean),
    depth = if ("depth_by_sequencing" %in% names(.)) {
      ifelse(
        !is.na(depth_by_sequencing) & depth_by_sequencing != "",
        as.character(depth_by_sequencing),
        str_extract(sample_base, "20M|50M|100M")
      )
    } else {
      str_extract(sample_base, "20M|50M|100M")
    },
    sample_label = vapply(as.character(sample_base), build_proportion_label, character(1)),
    sample_key   = ifelse(!is.na(depth), paste0(depth, "_", sample_label), paste0("noDepth_", sample_label))
  )

# Average predicted abundance per strain per sample per tool
mean_strain_abundance_no_K12_sakai <- master2_no_K12_sakai %>%
  group_by(
    strain_profiling_tool,
    sample_key,
    sample_label,
    ecoli_strain,
    phylogroup
  ) %>%
  summarise(
    mean_predicted_abundance = mean(predicted_abundance, na.rm = TRUE),
    .groups = "drop"
  )

# Phylogroup colour setup following previous set up (Dark2 palette)
phylo_levels <- c("A", "B1", "B2", "C", "D", "E", "F", "G")
highlight_phylo <- c("A", "E")
other_phylo_alpha <- 0.5

dark2 <- brewer.pal(8, "Dark2")
col_A_dark_exact <- dark2[2]
col_E_exact <- dark2[3]

lighten <- function(col, factor = 0.35) {
  rgbc <- col2rgb(col) / 255
  blended <- rgbc + (1 - rgbc) * factor
  rgb(t(blended), maxColorValue = 1)
}

col_A_light <- lighten(col_A_dark_exact, factor = 0.45)

make_shades <- function(base_color, n) {
  if (n == 1) return(base_color)
  colorRampPalette(c(base_color, lighten(base_color, 0.35)))(n)
}

remaining_phylo <- setdiff(phylo_levels, highlight_phylo)
remaining_colors <- dark2[setdiff(seq_along(dark2), c(2, 3))]

if (length(remaining_phylo) > length(remaining_colors)) {
  remaining_colors <- colorRampPalette(remaining_colors)(length(remaining_phylo))
} else {
  remaining_colors <- remaining_colors[seq_len(length(remaining_phylo))]
}
names(remaining_colors) <- remaining_phylo

# Labels, depth, and stacking order
mean_strain_abundance_no_K12_sakai <- mean_strain_abundance_no_K12_sakai %>%
  mutate(
    ecoli_label = paste0(ecoli_strain, " (", phylogroup, ")"),
    depth = factor(str_extract(sample_key, "20M|50M|100M"), levels = c("20M", "50M", "100M")),
    phylo_stack_order = case_when(
      phylogroup %in% setdiff(phylo_levels, highlight_phylo) ~ 1,
      phylogroup == "A" ~ 2,
      phylogroup == "E" ~ 3
    )
  )

stack_levels <- mean_strain_abundance_no_K12_sakai %>%
  distinct(ecoli_label, phylo_stack_order, phylogroup) %>%
  arrange(phylo_stack_order, factor(phylogroup, levels = phylo_levels), ecoli_label) %>%
  pull(ecoli_label)

mean_strain_abundance_no_K12_sakai$ecoli_label <-
  factor(mean_strain_abundance_no_K12_sakai$ecoli_label, levels = stack_levels)

# Build colour lookup by phylogroup
phylo_to_labels <- mean_strain_abundance_no_K12_sakai %>%
  distinct(ecoli_label, phylogroup) %>%
  arrange(factor(phylogroup, levels = phylo_levels)) %>%
  group_by(phylogroup) %>%
  summarise(labels = list(ecoli_label), .groups = "drop")

colour_lookup <- character(0)

for (i in seq_len(nrow(phylo_to_labels))) {
  phy <- phylo_to_labels$phylogroup[i]
  labs <- phylo_to_labels$labels[[i]]
  nlabs <- length(labs)
  
  if (phy == "A") {
    shades_final <- character(nlabs)
    mid_shades <- make_shades(col_A_dark_exact, max(1, nlabs))
    
    for (j in seq_along(labs)) {
      lab_j <- labs[j]
      if (grepl("REL606", lab_j, ignore.case = TRUE)) {
        shades_final[j] <- col_A_dark_exact
      } else if (grepl("\\bHS\\b", lab_j, ignore.case = TRUE) || grepl("HS_", lab_j, ignore.case = TRUE)) {
        shades_final[j] <- col_A_light
      } else {
        shades_final[j] <- mid_shades[((j - 1) %% length(mid_shades)) + 1]
      }
    }
    
  } else if (phy == "E") {
    shades_final <- if (nlabs == 1) col_E_exact else make_shades(col_E_exact, nlabs)
    
  } else {
    base_color <- remaining_colors[[as.character(phy)]]
    if (is.null(base_color)) base_color <- unname(sample(remaining_colors, 1))
    shades <- if (nlabs == 1) base_color else make_shades(base_color, nlabs)
    shades_final <- alpha(shades, other_phylo_alpha)
  }
  
  names(shades_final) <- labs
  colour_lookup <- c(colour_lookup, shades_final)
}

if (is.null(levels(mean_strain_abundance_no_K12_sakai$ecoli_label))) {
  mean_strain_abundance_no_K12_sakai <- mean_strain_abundance_no_K12_sakai %>%
    mutate(ecoli_label = factor(ecoli_label))
}

colour_lookup <- colour_lookup[levels(mean_strain_abundance_no_K12_sakai$ecoli_label)]

missing_labels <- setdiff(
  levels(mean_strain_abundance_no_K12_sakai$ecoli_label),
  names(colour_lookup)
)

if (length(missing_labels) > 0) {
  extra_cols <- colorRampPalette(dark2)(length(missing_labels))
  names(extra_cols) <- missing_labels
  colour_lookup <- c(colour_lookup, extra_cols)
  colour_lookup <- colour_lookup[levels(mean_strain_abundance_no_K12_sakai$ecoli_label)]
}

# Normalise within each facet group for stacked bars
plot_data <- mean_strain_abundance_no_K12_sakai %>%
  group_by(strain_profiling_tool, sample_label, depth) %>%
  mutate(
    total = sum(mean_predicted_abundance),
    mean_predicted_abundance = ifelse(total > 0, mean_predicted_abundance / total, 0)
  ) %>%
  select(-total) %>%
  ungroup()

# Plot
p_all_strains_no_K12_sakai_mod <-
  ggplot(plot_data,
         aes(x = sample_label, y = mean_predicted_abundance, fill = ecoli_label)) +
  geom_col(color = "black", width = 0.8, position = position_stack(reverse = FALSE)) +
  geom_text(
    aes(label = ifelse(mean_predicted_abundance > 0.08,
                       sprintf("%.2f", mean_predicted_abundance), "")),
    position = position_stack(vjust = 0.5),
    size = 2.4,
    color = "white",
    fontface = "bold"
  ) +
  facet_grid(
    rows = vars(depth),
    cols = vars(strain_profiling_tool),
    scales = "fixed",
    labeller = labeller(strain_profiling_tool = function(x) str_wrap(x, width = 16))
  ) +
  scale_fill_manual(
    values = colour_lookup,
    breaks = mean_strain_abundance_no_K12_sakai %>%
      distinct(ecoli_label, phylogroup) %>%
      arrange(factor(phylogroup, levels = phylo_levels), ecoli_label) %>%
      pull(ecoli_label),
    name = expression(italic("E. coli") ~ "strain (phylogroup)"),
      guide = guide_legend(ncol = 1)
    ) +
  labs(
    #title = "All predicted strains per sample (K12-MG1655 and O157:H7 str. Sakai excluded from reference database)",
    #subtitle = "Bars show mean predicted relative abundance across triplicates; legend shows strain (phylogroup)",
    x = "True K12-MG1655 : O157:H7 str. Sakai ratio in simulated metagenome sample",
    y = "Predicted relative abundance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 8),
    axis.title.x = element_text(vjust = -0.75),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key = element_rect(colour = "black", fill = NA),
    legend.text = element_text(size = 8)
  )

print(p_all_strains_no_K12_sakai_mod)
ggsave(plot = p_all_strains_no_K12_sakai_mod,
       filename = "all_strains_stacked_bars_no_K12_sakai.pdf",
       width = 14, height = 10)


#################################################################################################

### Fig. S4 - Specificity per tool when K12 and Sakai are removed from the reference database ###
# Counts/proportion of call types per strain-level profiler for samples ran without K12 or Sakai in the reference db

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(binom)

# Reload master table
master_no_K12_sakai <- read.csv("Simulated_metagenomes_K12_and_Sakai_removed_from_reference_database_raw_predicted_abundance.csv") # Available on FigShare

# Count per profiler × call category (include zeros)
count_table_no_K12_sakai <- master_no_K12_sakai %>%
  count(strain_profiling_tool, call_type, name = "n") %>%
  complete(strain_profiling_tool, call_type = c("TP","FP","FN","TN","other"), fill = list(n = 0))
# optional: put call categories in desired plotting order
#mutate(call_type = factor(call_type, levels = c("TP","FP","TN","FN","other")))

# Grouped bar plot (dodged bars: one bar per call type per profiler)
p_counts_no_K12_sakai <- ggplot(count_table_no_K12_sakai %>% filter(call_type %in% c("TP","FP","TN","FN")),
                                aes(x = strain_profiling_tool, y = n, fill = call_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, colour = "black", size = 0.2) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 3.5, fontface = "bold"
  ) +
  scale_y_continuous(name = "Count of calls", expand = expansion(mult = c(0,0.05))) +
  scale_fill_brewer(palette = "Dark2", name = "Call type") +
  labs(
    title = "Call-type counts by strain-level profiler (all samples; n = 8,316)",
    subtitle = "Counts of FP and TN across all samples; no threshold applied",
    x = "Strain-level profiling tool"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Print grouped bar plot - not used in FigS4 but informs it
print(p_counts_no_K12_sakai)
#ggsave("raw_calltype_counts_by_profiler_no_K12_sakai.png",
#      plot = p_counts_no_K12_sakai, width = 10, height = 6, dpi = 300)

# Continue to plot for specificity (TN and FP only; Fig. S4)

# Pivot to wide (TP, FP, TN, FN per profiler)
counts_wide_no_K12_sakai <- count_table_no_K12_sakai %>%
  filter(call_type %in% c("FP","TN")) %>%
  pivot_wider(names_from = call_type, values_from = n, values_fill = 0) %>%
  # ensure columns exist even if some profilers lack a call type
  mutate(
    FP = ifelse(is.na(FP), 0L, FP),
    TN = ifelse(is.na(TN), 0L, TN)
  )

# Compute specificity metric + 95% Wilson CI
metrics_table_no_K12_sakai <- counts_wide_no_K12_sakai %>%
  mutate(
    total_calls = FP + TN,
    specificity = ifelse(total_calls > 0, TN / total_calls, NA_real_)
  )

# Compute Wilson CIs
ci_table <- binom::binom.confint(
  x = metrics_table_no_K12_sakai$TN,
  n = metrics_table_no_K12_sakai$total_calls,
  methods = "wilson"
)

metrics_table_no_K12_sakai <- metrics_table_no_K12_sakai %>%
  mutate(
    ci_lower = ci_table$lower,
    ci_upper = ci_table$upper
  ) %>%
  arrange(desc(total_calls))

# Specificity figure (grouped bars): prepare long table for plotting
metrics_long_no_K12_sakai <- metrics_table_no_K12_sakai %>%
  select(strain_profiling_tool, specificity, ci_lower, ci_upper) %>%
  pivot_longer(cols = c(specificity),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("specificity"),
                         labels = c("Specificity")))

# position dodge object for consistent alignment of bars and labels
pd_no_K12_sakai <- position_dodge(width = 0.8)

# Plot
p_spec_no_K12_sakai <- ggplot(metrics_long_no_K12_sakai, aes(x = strain_profiling_tool, y = value, fill = strain_profiling_tool)) +
  geom_col(position = pd_no_K12_sakai, colour = "black", width = 0.7) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = pd_no_K12_sakai) +
  scale_y_continuous(name = "Specificity", limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Dark2", name = "") +
  labs(
    #title = "Specificity per strain-level profiler tool when K12-MG1655 and O157:H7 str. Sakai are excluded from the reference database",
    #subtitle = "Computed from raw FP and TN counts (n = 8,316); error bars show 95% Wilson confidence intervals",
    x = "Strain-level profiling tool",
    y = "Specificity"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  # numeric labels on top of each bar
  geom_text(
    data = metrics_long_no_K12_sakai,
    aes(label = ifelse(is.na(value), "NA", round(value, 3))),
    position = pd_no_K12_sakai,
    vjust = -1.5,
    size = 4,
    fontface = "bold"
  )

print(p_spec_no_K12_sakai)
ggsave("specificity_by_profiler_no_K12_sakai.pdf",
       plot = p_spec_no_K12_sakai, width = 12, height = 6)