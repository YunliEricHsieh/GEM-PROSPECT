setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

topDir <- "Results/"

# read & clean one table
read_and_clean <- function(path) {
  df <- read.table(path, header = TRUE, sep = ",")
  # Replace infinite with NA, negative with 0
  df[-1] <- lapply(df[-1], function(x) ifelse(is.infinite(x), NA, x))
  df[df < 0] <- 0
  # Drop rows where all numeric cols are NA
  df[rowSums(!is.na(df[-1])) > 0, ]
  # Drop cols where all numeric rows are NA
  df[colSums(!is.na(df[-1])) > 0, ]
}

# log10 fold-change analysis
apply_foldchange_filter <- function(df, sampling, threshold) {
  ids <- df[[1]]
  num_mat <- as.matrix(df[-1])
  mean_fluxes <- sampling$meanFlux[match(colnames(num_mat), sampling$RxnIndex)]
  mean_flux_mat <- matrix(mean_fluxes, nrow = nrow(num_mat), ncol = ncol(num_mat), byrow = TRUE)
  foldChange <- abs(num_mat) / abs(mean_flux_mat)
  
  # Identify cells to zero out: Not NA, and log10(foldChange) < threshold
  to_zero <- !is.na(foldChange) & (log10(foldChange) < threshold)
  num_mat[to_zero] <- 0
  result_df <- data.frame(ID = ids, num_mat)
  colnames(result_df)[1] <- colnames(df)[1]
  return(result_df)
}

# Phenotype data
phynotype_data <- read.table("Data/Mutant_phenotypes_table_filtered_final.csv", header = TRUE, sep = ",")

# Reversible and irreversible data
max_re <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re.csv"))
max_ir <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8.csv"))

# Sampling data
auto_re <- read.table(paste0(topDir, "flux_sampling/auto_sampling_re.csv"), header = TRUE, sep = ",")
auto_ir <- read.table(paste0(topDir, "flux_sampling/auto_sampling.csv"),    header = TRUE, sep = ",")

# Apply fold-change filter and combine
filtered_re <- apply_foldchange_filter(max_re, auto_re, -1)
filtered_ir <- apply_foldchange_filter(max_ir, auto_ir, -1)
T2_auto_max_data <- full_join(filtered_re, filtered_ir, by = "EnzymeID")

# Combine Raw max data for comparison
max_data <- full_join(max_re, max_ir, by = "EnzymeID")

phynotype_data <- phynotype_data %>% 
  filter(UniProtID %in% T2_auto_max_data$EnzymeID) %>%
  select(-GeneID, -mixo_NaCl_mixo, -mixo_P_mixo, -mixo_N_mixo)

# Rename columns
colnames(phynotype_data) <- c(
  "UniProtID", "auto/hetero", "auto/mixo", "auto/auto (CO2 3%)", 
  "mixo/mixo (CO2 3%)", "mixo/hetero", "mixo (Hypoosmotic 10%)/mixo", 
  "mixo (Hypoosmotic 25%)/mixo", "mixo (Hypoosmotic 75%)/mixo"
)

phynotype_data_melted <- phynotype_data %>%
  pivot_longer(
    cols = -UniProtID, 
    names_to = "Screens", 
    values_to = "Biomass_ratio"
  )

Zero_flux_count <- data.frame(
  Counts = log2(colSums(!is.na(max_data[, -1]) & (max_data[, -1] == 0))+1)
)

zero_flux_count_fold_change <- data.frame(
  Counts = log2(colSums(!is.na(T2_auto_max_data[, -1]) & (T2_auto_max_data[, -1] == 0))+1)
)

custom_theme <- theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"), # fill must be NA, not 'transparent'
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 17),
    plot.title = element_text(size = 18, face = "bold")
  )

plot_a <- ggplot(phynotype_data_melted, aes(x = Biomass_ratio, fill = Screens)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(
    title = "a.",
    x = "log2-transformed abundance ratio",
    y = "Density",
    fill = "Screens"
  ) +
  custom_theme +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

plot_b <- ggplot(Zero_flux_count, aes(y = Counts)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "b.",
    x = "Number of reactions with zero flux\nbefore fold-change analysis",
    y = "log2-transformed\nnumber of mutants"
  ) +
  custom_theme

plot_c <- ggplot(zero_flux_count_fold_change, aes(y = Counts)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(
    title = "c.",
    x = "Number of reactions with zero flux\nafter fold-change analysis",
    y = "log2-transformed\nnumber of mutants"
  ) +
  custom_theme

final_plot <- plot_a + (plot_b / plot_c) + 
  plot_layout(widths = c(2, 1))

print(final_plot)

ggsave("Results/figures/Fig 2.svg", final_plot, width = 17, height = 8)
