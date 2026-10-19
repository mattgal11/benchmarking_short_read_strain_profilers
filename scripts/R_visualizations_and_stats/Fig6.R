# ============================================================
# Fig. 6 Detection metrics of tools across 42 real-world 
# metagenomes with identical 10-strain spike-in
# ============================================================

# ============================================================
# Fig. 6a: Summary radar plot of detection metrics across 42 
# real-world metagenomes with identical 10-strain spike-in
# Detection metrics: sensitivity, specificity, precision, F1
# ============================================================

# ------------------------------------------------------------
# 1. Load libraries (ensure they are installed on your system)
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(fmsb)
library(RColorBrewer)
library(readxl)

# ------------------------------------------------------------
# 2. Load master final benchmarking dataset
# ------------------------------------------------------------
master <- read.csv("real_world_spike_ins_benchmarking_results.csv") # <-- Available on FigShare.

# ------------------------------------------------------------
# 2b. Normalise Strainify abundance estimates within each
# sample so that total predicted abundance = 1
# Based on reported solver issue
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

# Proof check normalisation for Strainify
print(master %>%
  filter(strain_profiling_tool == "Strainify") %>%
  group_by(sample_id) %>%
  summarise(
    total = sum(predicted_abundance, na.rm = TRUE)), n = 42)

# ------------------------------------------------------------
# 3. Prepare master dataset
# ------------------------------------------------------------
master2 <- master %>%
  mutate(
    predicted_abundance = as.numeric(predicted_abundance),
    # Assign each reference genome to one of the 10 true strains.
    # All other E. coli genomes are grouped into "Other".
    strain_category = case_when(
      strain_id %in% names(spiked_strains) ~
        unname(spiked_strains[strain_id]),
      TRUE ~ "Other"))

# ------------------------------------------------------------
# 4. Count TP / FP / TN / FN per profiling tool
# ------------------------------------------------------------
count_table <- master2 %>%
  count(
    strain_profiling_tool,
    call_type,
    name = "n"
  ) %>%
  complete(
    strain_profiling_tool,
    call_type = c("TP", "FP", "TN", "FN", "other"),
    fill = list(n = 0)
  ) %>%
  mutate(
    call_type = factor(
      call_type,
      levels = c("TP", "FP", "TN", "FN", "other")))

# ------------------------------------------------------------
# 5. Pivot to wide format
# ------------------------------------------------------------
counts_wide <- count_table %>%
  filter(
    call_type %in% c("TP", "FP", "TN", "FN")
  ) %>%
  pivot_wider(
    names_from = call_type,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    TP = ifelse(is.na(TP), 0L, TP),
    FP = ifelse(is.na(FP), 0L, FP),
    FN = ifelse(is.na(FN), 0L, FN),
    TN = ifelse(is.na(TN), 0L, TN))

# ------------------------------------------------------------
# 6. Calculate detection metrics
# ------------------------------------------------------------
metrics_table <- counts_wide %>%
  mutate(
    sensitivity = ifelse(
      (TP + FN) > 0,
      TP / (TP + FN),
      NA_real_),
    
    specificity = ifelse(
      (TN + FP) > 0,
      TN / (TN + FP),
      NA_real_),
    
    precision = ifelse(
      (TP + FP) > 0,
      TP / (TP + FP),
      NA_real_),
    
    F1 = ifelse(
      (TP + FP + FN) > 0,
      (2 * TP) / (2 * TP + FP + FN),
      NA_real_),
    
    total_calls = TP + FP + FN + TN
  ) %>%
  arrange(desc(total_calls))

print(metrics_table)

# ------------------------------------------------------------
# 7. Prepare radar plot data
# ------------------------------------------------------------
df_wide <- metrics_table %>%
  select(
    strain_profiling_tool,
    sensitivity,
    specificity,
    precision,
    F1
  ) %>%
  mutate(
    across(
      c(sensitivity, specificity, precision, F1),
      ~ ifelse(is.na(.), 0, as.numeric(.))))

# ------------------------------------------------------------
# 8. Build max / min rows
# ------------------------------------------------------------
max_row <- data.frame(
  sensitivity = 1,
  specificity = 1,
  precision = 1,
  F1 = 1)

min_row <- data.frame(
  sensitivity = 0,
  specificity = 0,
  precision = 0,
  F1 = 0)

radar_df <- rbind(
  max_row,
  min_row,
  df_wide %>%
    select(-strain_profiling_tool))

rownames(radar_df) <- c(
  "Max",
  "Min",
  df_wide$strain_profiling_tool)

# ------------------------------------------------------------
# 9. Colours (consistent with previous figures)
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

# ------------------------------------------------------------
# 10. Radar plot
# ------------------------------------------------------------
caxislabels <- c(
  "0",
  "0.25",
  "0.50",
  "0.75",
  "1.00")

pdf("Fig6a.pdf",
  width = 12,
  height = 6)

op <- par(
  mfrow = c(2, 4),
  mar = c(0.2, 3.5, 0.2, 2),
  oma = c(0, 0, 3, 0))


for (i in 3:nrow(radar_df)) {
  
  tool_name <- rownames(radar_df)[i]
  tool_colour <- tool_colours[tool_name]
  
  par(
    plt = c(0.10, 0.90, 0.15, 0.85))
  
  radarchart(
    radar_df[c(1, 2, i), ],
    axistype = 1,
    seg = 4,
    
    pcol = tool_colour,
    pfcol = scales::alpha(
      tool_colour,
      0.5),
    
    plwd = 1.8,
    plty = 1,
    
    cglcol = "grey90",
    cglty = 1,
    cglwd = 0.8,
    
    axislabcol = "grey40",
    
    vlcex = 1.0,
    
    caxislabels = caxislabels,
    
    vlabels = c(
      "sensitivity", "", "precision", "F1"),
    
    title = "")
  
  text(
    x = 0, y = 1.5,
    labels = tool_name,
    xpd = NA, cex = 1.2, font = 2)
  
  text(
    x = -1.05, y = 0,
    labels = "specificity",
    xpd = NA,
    adj = c(1, 0.5),
    cex = 1.0)
}

dev.off()
par(op)

###############################################################################

# ============================================================
# Fig. 6b: Detection threshold analysis. Assessing how 
# sensitivity, specificity, precision, and F1 score change 
# under detection thresholds (relative abundance cutoffs)
# ============================================================

# ------------------------------------------------------------
# 1. Detection-threshold analysis
# ------------------------------------------------------------
cutoffs <- c(
  0,
  0.0001,
  0.0005,
  0.001,
  0.005)

# ------------------------------------------------------------
# 2. Recalculate call type at each threshold
# ------------------------------------------------------------
master_all_threshold_cutoffs_list <- lapply(
  cutoffs,
  function(cutoff) {
    
    tmp <- master
    tmp$predicted_abundance <- as.numeric(
      tmp$predicted_abundance)
    tmp$predicted_abundance_cutoff <-
      tmp$predicted_abundance
    tmp$predicted_abundance_cutoff[
      tmp$predicted_abundance_cutoff < cutoff
    ] <- 0
    tmp$detected <-
      tmp$predicted_abundance_cutoff > 0
    tmp$call_type <- ifelse(
      tmp$strain_present_in_sample,
      ifelse(
        tmp$detected,
        "TP",
        "FN"),
      ifelse(
        tmp$detected,
        "FP",
        "TN"))
    
    tmp$threshold <- cutoff
    tmp
  }
)

combined_all_thresholds <-
  bind_rows(
    master_all_threshold_cutoffs_list)

# ------------------------------------------------------------
# 3. Count TP / FP / TN / FN at each threshold
# ------------------------------------------------------------
count_table_all_threshold_cutoffs <-
  combined_all_thresholds %>%
  
  mutate(
    call_type = as.character(call_type),
    call_type = trimws(call_type)
  ) %>%
  
  count(
    strain_profiling_tool,
    threshold,
    call_type,
    name = "n"
  ) %>%
  
  complete(
    strain_profiling_tool,
    threshold,
    call_type = c(
      "TP",
      "FP",
      "TN",
      "FN",
      "other"
    ),
    fill = list(n = 0)
  ) %>%
  
  mutate(
    call_type = factor(
      call_type,
      levels = c(
        "TP",
        "FP",
        "TN",
        "FN",
        "other"
      )
    )
  )

# ------------------------------------------------------------
# 4. Pivot to wide
# ------------------------------------------------------------
counts_wide_all_threshold_cutoffs <-
  count_table_all_threshold_cutoffs %>%
  
  filter(
    call_type %in% c(
      "TP",
      "FP",
      "TN",
      "FN"
    )
  ) %>%
  
  pivot_wider(
    names_from = call_type,
    values_from = n,
    values_fill = 0
  ) %>%
  
  mutate(
    TP = ifelse(is.na(TP), 0L, TP),
    FP = ifelse(is.na(FP), 0L, FP),
    FN = ifelse(is.na(FN), 0L, FN),
    TN = ifelse(is.na(TN), 0L, TN))

# ------------------------------------------------------------
# 4. Calculate detection metrics at each threshold
# ------------------------------------------------------------
set.seed(123)
n_boot <- 1000

metrics_table_all_threshold_cutoffs <-
  counts_wide_all_threshold_cutoffs %>%
  
  rowwise() %>%
  mutate(
    
    sensitivity = ifelse(
      (TP + FN) > 0,
      TP / (TP + FN),
      NA_real_),
    
    specificity = ifelse(
      (TN + FP) > 0,
      TN / (TN + FP),
      NA_real_),
    
    precision = ifelse(
      (TP + FP) > 0,
      TP / (TP + FP),
      NA_real_),
    
    F1 = ifelse(
      (TP + FP + FN) > 0,
      (2 * TP) /
        (2 * TP + FP + FN),
      NA_real_),
    
    # --------------------------------------------------------
    # Wilson 95% CI for sensitivity
    # --------------------------------------------------------
    sens_low = ifelse(
      (TP + FN) > 0,
      binom::binom.confint(
        TP,
        TP + FN,
        method = "wilson"
      )$lower,
      NA_real_),
    
    sens_high = ifelse(
      (TP + FN) > 0,
      binom::binom.confint(
        TP,
        TP + FN,
        method = "wilson"
      )$upper,
      NA_real_),
    
    # --------------------------------------------------------
    # Wilson 95% CI for specificity
    # --------------------------------------------------------
    spec_low = ifelse(
      (TN + FP) > 0,
      binom::binom.confint(
        TN,
        TN + FP,
        method = "wilson"
      )$lower,
      NA_real_),
    
    spec_high = ifelse(
      (TN + FP) > 0,
      binom::binom.confint(
        TN,
        TN + FP,
        method = "wilson"
      )$upper,
      NA_real_),
    
    # --------------------------------------------------------
    # Wilson 95% CI for precision
    # --------------------------------------------------------
    prec_low = ifelse(
      (TP + FP) > 0,
      binom::binom.confint(
        TP,
        TP + FP,
        method = "wilson"
      )$lower,
      NA_real_),
    
    prec_high = ifelse(
      (TP + FP) > 0,
      binom::binom.confint(
        TP,
        TP + FP,
        method = "wilson"
      )$upper,
      NA_real_)
    
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 4b. Bootstrap 95% CI for F1
# ------------------------------------------------------------
bootstrap_f1_ci <- function(
    tp, fp, fn,
    nboot = 1000) {
  
  if ((tp + fp + fn) == 0) {
    return(c(NA, NA))
  }
  
  counts <- c(
    rep("TP", tp),
    rep("FP", fp),
    rep("FN", fn)
  )
  
  boot_f1 <- replicate(
    nboot,
    {
      samp <- sample(
        counts,
        replace = TRUE
      )
      
      tp_b <- sum(samp == "TP")
      fp_b <- sum(samp == "FP")
      fn_b <- sum(samp == "FN")
      
      if (
        (2 * tp_b + fp_b + fn_b) == 0
      ) {
        return(NA)
      }
      
      (2 * tp_b) /
        (2 * tp_b + fp_b + fn_b)
    }
  )
  
  quantile(
    boot_f1,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
}


f1_ci_vals <- purrr::pmap(
  list(
    metrics_table_all_threshold_cutoffs$TP,
    metrics_table_all_threshold_cutoffs$FP,
    metrics_table_all_threshold_cutoffs$FN
  ),
  ~ bootstrap_f1_ci(
    ..1,
    ..2,
    ..3,
    nboot = n_boot
  )
)

metrics_table_all_threshold_cutoffs$F1_low <-
  purrr::map_dbl(
    f1_ci_vals,
    1)

metrics_table_all_threshold_cutoffs$F1_high <-
  purrr::map_dbl(
    f1_ci_vals,
    2)

print(metrics_table_all_threshold_cutoffs, n = 25)
# Save metrics_table_all_threshold_cutoffs to CSV for record
write.csv(metrics_table_all_threshold_cutoffs,
          file.path("real_world_spike_ins_detection_metrics_all_thresholds.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# 5. Long-format threshold data
# ------------------------------------------------------------
metrics_long_all_threshold_cutoffs <-
  metrics_table_all_threshold_cutoffs %>%
  
  select(
    strain_profiling_tool,
    threshold,
    sensitivity,
    specificity,
    precision,
    F1,
    sens_low,
    sens_high,
    spec_low,
    spec_high,
    prec_low,
    prec_high,
    F1_low,
    F1_high
  ) %>%
  
  pivot_longer(
    cols = c(
      sensitivity,
      specificity,
      precision,
      F1),
    
    names_to = "metric",
    values_to = "value"
  ) %>%
  
  mutate(
    ci_low = case_when(
      metric == "sensitivity" ~ sens_low,
      metric == "specificity" ~ spec_low,
      metric == "precision"   ~ prec_low,
      metric == "F1"          ~ F1_low),
    
    ci_high = case_when(
      metric == "sensitivity" ~ sens_high,
      metric == "specificity" ~ spec_high,
      metric == "precision"   ~ prec_high,
      metric == "F1"          ~ F1_high),
    
    metric = factor(
      metric,
      levels = c(
        "sensitivity",
        "specificity",
        "precision",
        "F1"
      ),
      
      labels = c(
        "Sensitivity",
        "Specificity",
        "Precision",
        "F1 score"
      )
    )
  )

# ------------------------------------------------------------
# 6. Format threshold labels
# ------------------------------------------------------------
metrics_long_all_threshold_cutoffsv2 <-
  metrics_long_all_threshold_cutoffs %>%
  
  mutate(
    threshold_f = factor(
      threshold,
      levels = cutoffs,
      labels = {
        lbl <- sub(
          "\\.?0+$",
          "",
          format(
            cutoffs,
            scientific = FALSE
          )
        )
        lbl[lbl == ""] <- "0"
        lbl
      }
    )
  )

# ------------------------------------------------------------
# 7. Define Dark2 colours (PanTax absent) and plot detection 
# metrics across thresholds
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

p_final_benchmark_threshold_cutoffs <-
  ggplot(metrics_long_all_threshold_cutoffsv2,
    aes(x = threshold_f, y = value,
      color = strain_profiling_tool,
      group = strain_profiling_tool)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15, alpha = 0.6) +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))) +
  scale_colour_manual(values = tool_colours) +
  labs(
    x = "Detection threshold (relative abundance cutoff)",
    y = "Detection metric score",
    color = "Strain-level profiling tool") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 14)),
    axis.text.x = element_text(angle = 37.5, vjust = 0.5, size = 10, face = "bold"))

print(p_final_benchmark_threshold_cutoffs)
ggsave("Fig6b.pdf", p_final_benchmark_threshold_cutoffs, width = 14, height = 9)

# ============================================================
# Cowplot combined Fig. 6 for final publication
# ============================================================

# (Re)load libraries
library(cowplot)
library(ggplotify)
library(fmsb)
library(scales)
library(ggplot2)

# Need to wrap the radar plot as a ggplot object first otherwise cowplot will not work
radar_ggplot <- as.ggplot(~{
#op <- par(
#  mfrow = c(4, 2),
#  mar = c(0.2, 3.5, 0.2, 2),
#  oma = c(0, 0, 3, 0))

op <- par(
  mfrow = c(4, 2),
  mar = c(0, 0.5, 0, 0.5),
  oma = c(0, 0, 0, 0))

for (i in 3:nrow(radar_df)) {
  
  tool_name <- rownames(radar_df)[i]
  tool_colour <- tool_colours[tool_name]
  
  par(
    plt = c(0.05, 0.95, 0.10, 0.90))
  
  radarchart(
    radar_df[c(1, 2, i), ],
    axistype = 1,
    seg = 4,
    
    pcol = tool_colour,
    pfcol = scales::alpha(
      tool_colour,
      0.5),
    
    plwd = 1.8,
    plty = 1,
    
    cglcol = "grey90",
    cglty = 1,
    cglwd = 0.8,
    
    axislabcol = "grey40",
    
    vlcex = 1.0,
    
    caxislabels = caxislabels,
    
    vlabels = c(
      "sensitivity", "", "precision", "F1"),
    
    title = "")
  
  text(
    x = 0, y = 1.5,
    labels = tool_name,
    xpd = NA, cex = 1.2, font = 2)
  
  text(
    x = -1.05, y = 0,
    labels = "specificity",
    xpd = NA,
    adj = c(1, 0.5),
    cex = 1.0)
    }
  }
)

print(radar_ggplot)

# Fig.6 combined with cowplot
p6a <- radar_ggplot
p6b <- p_final_benchmark_threshold_cutoffs

plot_grid(p6a, p6b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(1.3, 1.7),
          align = "h", axis = "tb")
ggsave("combined_Fig6.pdf", width = 12, height = 8)