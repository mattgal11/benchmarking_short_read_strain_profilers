# This script produces the three figures from the paper associated with the ZymoBIOMICS® D6331 gut microbiome standard dataset

# load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# load in the data and remove PathoScope final guess rows (see Methods)
zymo_error <- read.csv("Zymobiomics_D6331_raw_predicted_abundance.csv") # available on FigShare
zymo_error <- zymo_error[-c(16:20), ]

# data prep
zymo_error$Predicted_abundance <- as.numeric(zymo_error$Predicted_abundance)
zymo_error$Absolute_error <- as.numeric(zymo_error$Absolute_error)
zymo_error$Error <- as.numeric(zymo_error$Error)

# set the required order of tools for consistency across figures
tool_order <- c(
  "Ground truth",
  "PanTax",
  "PathoScope",
  "StrainGE",
  "Strainify",
  "StrainR2",
  "StrainScan"
)

# set plotting colours of Dark2 to be carried through all figures
# 6 required for this dataset - these will be the colours of each tool
tool_cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02")

# set plotting colours for strains in the ZYMO dataset
strain_cols <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7")

# --- Fig. 1a. Predicted relative abundance shown as stacked bars ---

# make a ground-truth bar with the same strain proportions
gt_bar <- zymo_error %>%
  distinct(Ecoli_strain) %>%
  mutate(
    Strain_profiling_tool = "Ground truth",
    Predicted_abundance = 0.20
  )

# combine with predictions
zymo_error_stack <- bind_rows(zymo_error, gt_bar) %>%
  mutate(
    Strain_profiling_tool = factor(Strain_profiling_tool, levels = tool_order)
  )

p_zymo_pred_abundance_stacked <- ggplot(
  zymo_error_stack,
  aes(x = Strain_profiling_tool, y = Predicted_abundance, fill = Ecoli_strain)
  ) +
  geom_col(width = 0.8, colour = NA) +
  geom_text(
    aes(label = sprintf("%.2f", Predicted_abundance)),
    position = position_stack(vjust = 0.5),
    size = 4, angle = 0, fontface = "bold", show.legend = FALSE, color = "white"
  ) +
  scale_fill_manual(values = strain_cols, name = expression(italic(E.~coli) ~ "strain")) +
  labs(
    x = "Strain-level profiling tool",
    y = "Predicted relative abundance"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 0, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 15)),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(t = 10, r = 50, b = 30, l = 10)
  )

# print and save
print(p_zymo_pred_abundance_stacked)
ggsave("Zymo_known_comparisons_figures/ZYMO_predicted_abundance_stacked_bars_updated_version.pdf",
       plot = p_zymo_pred_abundance_stacked, width = 14, height = 10)


# --- Fig. 1b. Mean absolute error per tool ---

#Plot absolute error as a mean per tool to show overall accuracy

# summarise mean + SD per tool
zymo_absolute_error_summary <- zymo_error %>%
  group_by(Strain_profiling_tool) %>%
  summarise(
    MAE = mean(Absolute_error, na.rm = TRUE),
    SD  = sd(Absolute_error, na.rm = TRUE),
    n   = n(),
    .groups = "drop")

# plotting
p_ZYMO_mean_absolute_error <- ggplot(zymo_absolute_error_summary, aes(x = Strain_profiling_tool, y = MAE, fill = Strain_profiling_tool)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(aes(ymin = MAE - SD, ymax = MAE + SD), width = 0.2, linewidth = 0.6) +
  geom_text(aes(y = MAE + SD, label = round(MAE, 3), vjust = -0.8),
          size = 3.5, fontface = "bold") +
  scale_fill_manual(values = tool_cols) +
  labs(x = "Strain-level profiling tool",
      y = "Mean absolute error"
      # OPTIONAL graph labels (unhash)
    #title = "Mean absolute error per strain-level profiling tool",
    #subtitle = "Error bars show standard deviation from predicted abundance guesses across E. coli strains (n=5)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.2, size = 9),
        axis.title.x = element_text(margin = margin(t = 13))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# print and save
print(p_ZYMO_mean_absolute_error)
ggsave("Zymo_known_comparisons_figures/ZYMO_absolute_error_per_tool_bar_plot.pdf",
       plot = p_ZYMO_mean_absolute_error, width = 12, height = 6)


### cowplot joint multi-panel Fig.1 for final publication ###

library(cowplot)

# Fig.1
p1a <- p_zymo_pred_abundance_stacked
p1b <- p_ZYMO_mean_absolute_error

plot_grid(p1a, p1b, labels = c("(a)", "(b)"), ncol = 2, rel_widths = c(1.75, 1.25))
ggsave("Zymo_known_comparisons_figures/ZYMO_Fig1_combined.pdf", width = 16, height = 8)


##################################################################################################


# --- Fig. S1. Error values from predicted abundance (included within this results narrative) ---

# aggregate if duplicates exist
zymo_plot <- zymo_error %>%
  group_by(Ecoli_strain, Strain_profiling_tool) %>%
  summarise(Error = mean(Error, na.rm = TRUE), .groups = "drop")

# ensure factor order for strains
zymo_plot$Ecoli_strain <- factor(zymo_plot$Ecoli_strain, levels = unique(zymo_plot$Ecoli_strain))

# plotting
dodge <- position_dodge(width = 0.9)
p_zymo_error_strain_facet <- ggplot(zymo_plot, aes(x = Strain_profiling_tool, y = Error, fill = Strain_profiling_tool)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = round(Error, 3),
                hjust = ifelse(Error >= 0, -0.3, 1.3)),
            position = position_dodge(width = 0.9),
            size = 2.5, fontface = "bold", angle = 90) +
  facet_wrap(~ Ecoli_strain, ncol = 6, scales = "fixed") +
  scale_fill_manual(values = tool_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  labs(x = "Strain-level profiling tool", y = "Error"
       # OPTIONAL graph titles (unhash)
       #title = "Error of tools in predicting relative abundance of E. coli strains in ZymoBIOMICS® D6331 Gut Microbiome Standard",
       #subtitle = "Faceted by strain to show under/overestimation of each tool"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 52.5, hjust = 1, size = 9),
        axis.title.x = element_text(margin = margin(t = 10)),
        strip.text = element_text(size = 10, face = "bold"))

# print and save ZYMO_error_bar_plot_strain_facet
print(p_zymo_error_strain_facet)
ggsave("Zymo_known_comparisons_figures/ZYMO_FigS1.pdf",
       plot = p_zymo_error_strain_facet, width = 12, height = 6)