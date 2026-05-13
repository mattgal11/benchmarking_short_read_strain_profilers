# Final plotting scripts for Fig. 3 (Detection metrics radar plot and thresholds)

### Fig. 3a - radar plot of sensitivity, specificity, precision and F1 score by strain-level profiler ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(fmsb)
library(RColorBrewer)

# 1) Count per profiler × call category (include zeros) #
master <- read.csv("Simulated_metagenomes_all_references_raw_predicted_abundances.csv") # dataset available on FigShare.
count_table <- master %>%
  count(strain_profiling_tool, call_type, name = "n") %>%
  # ensure all combinations exist (so missing categories show as zero)
  complete(strain_profiling_tool, call_type = c("TP","FP","FN","TN","other"), fill = list(n = 0)) %>%
  # optional: put call categories in desired plotting order
  mutate(call_type = factor(call_type, levels = c("TP","FP","TN","FN","other")))

# Optional: save raw counts table
 #write.csv(count_table, file.path("final_raw_call-type_counts_by_profiler.csv"), row.names = FALSE)

# 2) Pivot to wide (TP, FP, TN, FN per profiler; allows sensitivity/specificity/precision) #
counts_wide <- count_table %>%
  filter(call_type %in% c("TP","FP","TN","FN")) %>%
  pivot_wider(names_from = call_type, values_from = n, values_fill = 0) %>%
  # ensure columns exist even if some profilers lack a call type
  mutate(
    TP = ifelse(is.na(TP), 0L, TP),
    FP = ifelse(is.na(FP), 0L, FP),
    FN = ifelse(is.na(FN), 0L, FN),
    TN = ifelse(is.na(TN), 0L, TN)
  )

# 3) Compute metrics: sensitivity, specificity, precision, F1 #
metrics_table <- counts_wide %>%
  mutate(
    sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
    specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
    precision   = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
    F1          = ifelse((TP + FP + FN) > 0, (2 * TP) / (2 * TP + FP + FN), NA_real_),
    total_calls = TP + FP + FN + TN
  ) %>%
  arrange(desc(total_calls))

# Optional: write metrics to CSV 
 #write.csv(metrics_table, file.path("final_profiler_detection_metrics_per_tool.csv"), row.names = FALSE)

# 4) Plot radar chart of metrics (sensitivity, specificity, precision, F1) per profiler #

# Prepare wide table (ensure numeric, NAs -> 0)
df_wide <- metrics_table %>%
  select(strain_profiling_tool, sensitivity, specificity, precision, F1) %>%
  mutate(across(c(sensitivity, specificity, precision, F1),
                ~ ifelse(is.na(.), 0, as.numeric(.))))

# Build max/min rows as DATA FRAMES (0-1 scale)
max_row <- data.frame(sensitivity = 1, specificity = 1, precision = 1, F1 = 1)
min_row <- data.frame(sensitivity = 0, specificity = 0, precision = 0, F1 = 0)

radar_df <- rbind(max_row, min_row, df_wide %>% select(-strain_profiling_tool))
rownames(radar_df) <- c("Max", "Min", df_wide$strain_profiling_tool)

# Colour palette: use Dark2
n_tools <- nrow(df_wide)
cols <- brewer.pal(8, "Dark2")
cols_use <- rep(cols, length.out = n_tools)

# Axis label strings (0-1)
caxislabels <- c("0", "0.25", "0.50", "0.75", "1.00")

pdf("summary_radar_plot.pdf", width = 12, height = 6)

# set plotting layout before drawing
op <- par(
  mfrow = c(2,4),
  mar = c(0.2, 3.5, 0.2, 2),
  oma = c(0,0,3,0)
)

# Plot each profiler as a separate radar (Max, Min, and data row)
for (i in 3:nrow(radar_df)) {
  par(plt = c(0.10, 0.90, 0.15, 0.85))
  radarchart(
    radar_df[c(1, 2, i), ],
    axistype = 1,
    seg = 4,
    pcol = cols_use[i - 2],
    pfcol = alpha(cols_use[i - 2], 0.5),
    plwd = 1.8,
    plty = 1,
    cglcol = "grey90",
    cglty = 1,
    cglwd = 0.8,
    axislabcol = "grey40",
    vlcex = 1.0,
    caxislabels = caxislabels,
    vlabels = c("sensitivity", "", "precision", "F1"),
    title = ""
  )

  text(
    x = 0, y = 1.5,
    labels = rownames(radar_df)[i],
    xpd = NA,
    cex = 1.2,
    font = 2
  )
  
  text(
    x = -1.05, y = 0,
    labels = "specificity", # add in text label for specificity here to avoid overlap
    xpd = NA,
    adj = c(1, 0.5),
    cex = 1.0
  )
}

# Add global title and subtitle in the outer margin (OMA)
#mtext("Radar plot of sensitivity, specificity, precision and F1 score by strain-level profiler",
#      outer = TRUE, cex = 1.4, adj = 0.5, line = 2.3)
#mtext("Computed from raw TP, FP, TN and FN counts; metrics scaled 0-1",
#      outer = TRUE, cex = 1.0, adj = 0.5, line = 0.8)

dev.off()
par(op)  # restore old par settings if needed


### Fig. 3b - thresholds for detection metrics (sensitivity, specificity, precision, F1) per profiler ###

# (re)load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(binom)
library(purrr)

# 1) Calculate updated detection threshold .csv files (if needed - if not use those provided)

# parameters: change cutoffs here
cutoffs <- c(0, 0.0001, 0.0005, 0.001, 0.005)
out_dir <- "differential_abundance_threshold_cutoffs"
if (!dir.exists(out_dir)) dir.create(out_dir)

# reload data (if needed) and ensure required cols exist
master <- read.csv("Simulated_metagenomes_all_references_raw_predicted_abundances.csv") # available on FigShare.
required_cols <- c("strain_present_in_sample", "predicted_abundance", "call_type")
stopifnot(all(required_cols %in% colnames(master)))

# run loop for establishing different threshold cutoffs
for (c in cutoffs) {
  tmp <- master 
  
  ## apply cutoff: zero-out values below the cutoff (keep zeros)
  tmp$predicted_abundance_cutoff <- tmp$predicted_abundance
  tmp$predicted_abundance_cutoff[tmp$predicted_abundance_cutoff < c] <- 0
  
  ## define detected (logical) based on cutoff
  # use >= if you want to consider exact equals as detected; with zeroing >0 works equally
  tmp$detected <- tmp$predicted_abundance_cutoff > 0
  
  ## recalculate call_type using your Excel logic:
  # IF(strain_present_in_sample, IF(detected, "TP","FN"), IF(detected, "FP","TN"))
  tmp$call_type <- ifelse(tmp$strain_present_in_sample,
                          ifelse(tmp$detected, "TP", "FN"),
                          ifelse(tmp$detected, "FP", "TN"))
  
  ## save object (named) and CSV for record
  label <- gsub("\\.", "_", formatC(c, format = "f", digits = 6)) # e.g. 0.000500 -> "0_000500"
  varname <- make.names(paste0("master_cutoff_", label))
  assign(varname, tmp, envir = .GlobalEnv)
  
  fname <- file.path(out_dir, paste0(varname, ".csv"))
  write.csv(tmp, fname, row.names = FALSE)
  
  message("Saved object: ", varname, "  and CSV: ", fname)
}

# 2) Put master_cutoff data.frames and thresholds here

master_cutoff_0_000000 <- read.csv("Differential_abundance_threshold_cutoffs/master_cutoff_0_000000.csv")
master_cutoff_0_000100 <- read.csv("Differential_abundance_threshold_cutoffs/master_cutoff_0_000100.csv")
master_cutoff_0_000500 <- read.csv("Differential_abundance_threshold_cutoffs/master_cutoff_0_000500.csv")
master_cutoff_0_001000 <- read.csv("Differential_abundance_threshold_cutoffs/master_cutoff_0_001000.csv")
master_cutoff_0_005000 <- read.csv("Differential_abundance_threshold_cutoffs/master_cutoff_0_005000.csv")

thresholds <- c(0, 0.0001, 0.0005, 0.001, 0.005)

master_all_threshold_cutoffs_list <- list(
  master_cutoff_0_000000,  # threshold = 0
  master_cutoff_0_000100,  # threshold = 0.0001
  master_cutoff_0_000500,  # threshold = 0.0005
  master_cutoff_0_001000,  # threshold = 0.001
  master_cutoff_0_005000   # threshold = 0.005
)

# check lengths match
if(length(master_all_threshold_cutoffs_list) != length(thresholds)) {
  stop("length(master_all_threshold_cutoffs_list) must equal length(thresholds). Update names/vectors accordingly.")
}

# 3) Add a threshold column to each df and bind rows
# Use Map (safer/clearer than mapply with name shadowing). Ensure thresholds are numeric.
combined_all_thresholds <- Map(
  function(df, thr) {
    df %>% mutate(threshold = as.numeric(thr))
  },
  master_all_threshold_cutoffs_list,
  thresholds
) %>% bind_rows()

# 4) Count TP/FP/TN/FN per tool × threshold (include zeros)
# Ensure call_type is character (and trimmed) so matching works.
count_table_all_threshold_cutoffs <- combined_all_thresholds %>%
  mutate(call_type = as.character(call_type),
         call_type = trimws(call_type)) %>%
  count(strain_profiling_tool, threshold, call_type, name = "n") %>%
  complete(
    strain_profiling_tool,
    threshold,
    call_type = c("TP","FP","TN","FN","other"),
    fill = list(n = 0)
  ) %>%
  mutate(call_type = factor(call_type, levels = c("TP","FP","TN","FN","other")))

# 5) Pivot to wide (TP, FP, TN, FN per tool × threshold)
counts_wide_all_threshold_cutoffs <- count_table_all_threshold_cutoffs %>%
  filter(call_type %in% c("TP","FP","TN","FN")) %>%
  pivot_wider(names_from = call_type, values_from = n, values_fill = 0) %>%
  mutate(
    TP = ifelse(is.na(TP), 0L, TP),
    FP = ifelse(is.na(FP), 0L, FP),
    FN = ifelse(is.na(FN), 0L, FN),
    TN = ifelse(is.na(TN), 0L, TN)
  )

# 6) Compute metrics per tool × threshold

set.seed(123)  # reproducibility for bootstrap
n_boot <- 1000 # increase to 2000+ for final thesis figures

  # 6.1) Compute metrics + CIs

  metrics_table_all_threshold_cutoffs <- counts_wide_all_threshold_cutoffs %>%
  rowwise() %>%
  mutate(
    # Core metrics
    sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),
    specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),
    precision   = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),
    F1          = ifelse((TP + FP + FN) > 0,
                         (2 * TP) / (2 * TP + FP + FN),
                         NA_real_),
    
    # --- Wilson CIs (binomial) ---
    sens_low  = ifelse((TP + FN) > 0,
                       binom.confint(TP, TP + FN, method = "wilson")$lower,
                       NA_real_),
    sens_high = ifelse((TP + FN) > 0,
                       binom.confint(TP, TP + FN, method = "wilson")$upper,
                       NA_real_),
    
    spec_low  = ifelse((TN + FP) > 0,
                       binom.confint(TN, TN + FP, method = "wilson")$lower,
                       NA_real_),
    spec_high = ifelse((TN + FP) > 0,
                       binom.confint(TN, TN + FP, method = "wilson")$upper,
                       NA_real_),
    
    prec_low  = ifelse((TP + FP) > 0,
                       binom.confint(TP, TP + FP, method = "wilson")$lower,
                       NA_real_),
    prec_high = ifelse((TP + FP) > 0,
                       binom.confint(TP, TP + FP, method = "wilson")$upper,
                       NA_real_)
  ) %>%
  ungroup()

  # 6.2) Bootstrap CIs for F1 score

  bootstrap_f1_ci <- function(tp, fp, fn, nboot = 1000) {
  if ((tp + fp + fn) == 0) return(c(NA, NA))
  
  counts <- c(rep("TP", tp),
              rep("FP", fp),
              rep("FN", fn))
  
  boot_f1 <- replicate(nboot, {
    samp <- sample(counts, replace = TRUE)
    tp_b <- sum(samp == "TP")
    fp_b <- sum(samp == "FP")
    fn_b <- sum(samp == "FN")
    
    if ((2 * tp_b + fp_b + fn_b) == 0) return(NA)
    (2 * tp_b) / (2 * tp_b + fp_b + fn_b)
  })
  
  quantile(boot_f1, probs = c(0.025, 0.975), na.rm = TRUE)
  }

  f1_ci_vals <- pmap(
  list(metrics_table_all_threshold_cutoffs$TP,
       metrics_table_all_threshold_cutoffs$FP,
       metrics_table_all_threshold_cutoffs$FN),
  ~ bootstrap_f1_ci(..1, ..2, ..3, nboot = n_boot)
  )

metrics_table_all_threshold_cutoffs$F1_low  <- map_dbl(f1_ci_vals, 1)
metrics_table_all_threshold_cutoffs$F1_high <- map_dbl(f1_ci_vals, 2)

  # save metrics_table_all_threshold_cutoffs to CSV for record
  # write.csv(metrics_table_all_threshold_cutoffs,
  #          file.path("Differential_abundance_figures/near_final_pdf_versions/final_metrics_with_CIs_all_threshold_cutoffs.csv"),
  #          row.names = FALSE)

  # 6.3) Pivot to long with CIs

  metrics_long_all_threshold_cutoffs <- metrics_table_all_threshold_cutoffs %>%
  select(strain_profiling_tool, threshold,
         sensitivity, specificity, precision, F1,
         sens_low, sens_high,
         spec_low, spec_high,
         prec_low, prec_high,
         F1_low, F1_high) %>%
  pivot_longer(
    cols = c(sensitivity, specificity, precision, F1),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    ci_low = case_when(
      metric == "sensitivity" ~ sens_low,
      metric == "specificity" ~ spec_low,
      metric == "precision"   ~ prec_low,
      metric == "F1"          ~ F1_low
    ),
    ci_high = case_when(
      metric == "sensitivity" ~ sens_high,
      metric == "specificity" ~ spec_high,
      metric == "precision"   ~ prec_high,
      metric == "F1"          ~ F1_high
    ),
    metric = factor(metric,
                    levels = c("sensitivity","specificity","precision","F1"),
                    labels = c("Sensitivity","Specificity","Precision","F1 score"))
  )

# 7). Plotting script for detection threshold cutoffs

metrics_long_all_threshold_cutoffsv2 <- metrics_long_all_threshold_cutoffs %>%
  mutate(
    threshold_f = factor(
      threshold,
      levels = thresholds,
      labels = {
        lbl <- sub("\\.?0+$", "", format(thresholds, scientific = FALSE))
        lbl[lbl == ""] <- "0"
        lbl
      })
  )

p_all_metrics_all_threshold_cutoffsv2 <- ggplot(
  metrics_long_all_threshold_cutoffsv2,
  aes(x = threshold_f, y = value,
      color = strain_profiling_tool,
      group = strain_profiling_tool)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15,
                alpha = 0.6) +
  facet_wrap(~ metric, ncol = 2) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_color_brewer(palette = "Dark2",
                     name = "Strain-level profiling tool") +
  labs(
    #title = "Strain-level profiling tool performance across increasing detection thresholds",
    #subtitle = "Error bars represent 95% CIs (Wilson intervals for sensitivity, specificity and precision; bootstrap intervals for F1 score).",
    x = "Detection threshold (relative abundance cutoff)",
    y = "Detection metric score",
    color = "Strain-level profiling tool"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 14)),
    axis.text.x = element_text(angle = 37.5, vjust = 0.5, size = 10, face = "bold")
  )

print(p_all_metrics_all_threshold_cutoffsv2)

ggsave("line_plot_all_profiler_metrics_discrete_x_axis.pdf",
  p_all_metrics_all_threshold_cutoffsv2, width = 12, height = 6)


### cowplot combined figure 3 ###

# need to wrap the radar plot as a ggplot object first

library(cowplot)
library(ggplotify)
library(fmsb)
library(scales)
library(ggplot2)

radar_ggplot <- as.ggplot(~{
  op <- par(
    mfrow = c(4, 2),
    mar = c(0.8, 0.8, 0.8, 0.8),
    oma = c(0, 0, 0, 0),
    pty = "s"
  )
  par(xaxs = "i", yaxs = "i")
  on.exit(par(op), add = TRUE)
  
  for (i in 3:nrow(radar_df)) {
    radarchart(
      radar_df[c(1, 2, i), ],
      axistype = 1,
      seg = 4,
      pcol = cols_use[i - 2],
      pfcol = alpha(cols_use[i - 2], 0.5),
      plwd = 1.8,
      plty = 1,
      cglcol = "grey90",
      cglty = 1,
      cglwd = 0.8,
      axislabcol = "grey40",
      vlcex = 0.9,
      caxislabels = caxislabels,
      vlabels = c("sensitivity", "", "precision", "F1"),
      title = ""
    )
    
    text(
      x = -2, y = 0.85,
      labels = paste(strwrap(rownames(radar_df)[i], width = 15), collapse = "\n"),
      xpd = NA,
      adj = c(0, 0.5),
      cex = 1.0,
      font = 2
    )
    
    text(
      x = -1.08, y = 0,
      labels = "specificity",
      xpd = NA,
      adj = c(1, 0.5),
      cex = 0.85
    )
  }
})
print(radar_ggplot)

# Fig.3 combined with cowplot
p3a <- radar_ggplot + theme(aspect.ratio = 1)
p3b <- p_all_metrics_all_threshold_cutoffsv2

plot_grid(p3a, p3b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(1.4, 1.6),
  align = "h", axis = "tb")
ggsave("Simulated_metagenomes_Fig3.pdf", width = 12, height = 8)