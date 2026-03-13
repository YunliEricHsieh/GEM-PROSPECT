setwd('/Users/yunli/GEM-PROSPECT/')

phynotype_data <- read.table("Data/Mutant_phenotypes_table_filtered_final.csv", header = TRUE, sep = ",")

# read & clean one table
read_and_clean <- function(path) {
  df <- read.table(path, header=TRUE, sep=",")
  # infinite → NA, negative → 0
  df[-1] <- lapply(df[-1], function(x) ifelse(is.infinite(x), NA, x))
  df[df < 0] <- 0
  # drop rows where all numeric cols are NA
  df[rowSums(!is.na(df[-1])) > 0, ]
}

topDir <- "Results/"

# read reversible and irreversible data
max_re <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re.csv"))
max_ir <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8.csv"))

auto_re   = read.table(paste0(topDir, "flux_sampling/auto_sampling_re.csv"),   header=TRUE, sep=",")
auto_ir   = read.table(paste0(topDir, "flux_sampling/auto_sampling.csv"),   header=TRUE, sep=",")

# log10 fold-change analysis
apply_foldchange_filter <- function(df, sampling, threshold) {
  for (i in 2:ncol(df)){
    colid <- colnames(df)[i]
    rowindex <- which(sampling$RxnIndex == colid)
    meanValue <- sampling$meanFlux[rowindex]
    for (j in 1:nrow(df)){
      fluxValue <- df[j,i]
      if (!is.na(fluxValue)){
        foldChange <- abs(fluxValue) / abs(meanValue)
        if (foldChange != 0 && log10(foldChange) < threshold){
          df[j,i] <- 0
        }
      }
    }
  }
  return(df)
}

T2_auto_max_data <- cbind(apply_foldchange_filter(max_re, auto_re, -2),
                          apply_foldchange_filter(max_ir, auto_ir, -2)[-1])

library(dplyr)
phynotype_data <- phynotype_data %>% filter(UniProtID %in% T2_auto_max_data$EnzymeID)

phynotype_data <- phynotype_data %>% select(-GeneID, -mixo_NaCl_mixo, -mixo_P_mixo, -mixo_N_mixo)
colnames(phynotype_data) <- c("UniProtID", "auto/hetero", "auto/mixo", "auto/auto (CO2 3%)", "mixo/mixo (CO2 3%)", "mixo/hetero",
                              "mixo (Hypoosmotic 10%)/mixo", "mixo (Hypoosmotic 25%)/mixo", "mixo (Hypoosmotic 75%)/mixo")

library(ggplot2)

# melt the data for easier plotting
library(reshape2)
phynotype_data_melted <- melt(phynotype_data, id.vars = "UniProtID")
colnames(phynotype_data_melted) <- c("UniProtID", "Screens", "Biomass_ratio")

# plot the distribution for each column
biomass_density <- ggplot(phynotype_data_melted, aes(x = Biomass_ratio, fill = Screens)) +
    geom_density(alpha = 0.5) +
    scale_fill_brewer(palette = "Set3") +
    labs(title = "a.",
             x = "log2-transformed abundance ratio",
             y = "Density",
             fill = "Screens") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 18),
                plot.title = element_text(size = 18, face = "bold"),
                legend.text = element_text(size = 16),
                legend.title = element_text(size = 16))
print(biomass_density)

# number of mutant with 0 flux of rxns
max_data <- cbind(max_re, max_ir[-1])

Zero_flux_count <- data.frame(
  Counts = log2(colSums((!is.na(max_data[, -1]) & 
                      (max_data[, -1] == 0) ))),
  stringsAsFactors = FALSE
)

# Plot the distribution of zero flux count
zero_flux_distribution <- ggplot(Zero_flux_count, aes(y = Counts)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "b.",
         x = "Number of reactions with zero flux\nbefore fold-change analysis",
         y = "log2-transformed\nnumber of mutants") +
    theme_minimal() +
    theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
          axis.text = element_text(size = 17),
          axis.title.x = element_text(hjust = 0.5),
          axis.title = element_text(size = 17),
          plot.title = element_text(size = 17, face = "bold"))

print(zero_flux_distribution)

# number of mutant with 0 flux of rxns after fold-change analysis
zero_flux_count_fold_change <- data.frame(
  Counts = log2(colSums((!is.na(T2_auto_max_data[, -1]) & 
                      (T2_auto_max_data[, -1] == 0)))),
  stringsAsFactors = FALSE
)

# Plot the distribution of zero flux count
zero_flux_distribution_fold_change <- ggplot(zero_flux_count_fold_change, aes(y = Counts)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "c.",
       x = "Number of reactions with zero flux\nafter fold-change analysis",
       y = "log2-transformed\nnumber of mutants") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
        axis.text = element_text(size = 17),
        axis.title.x = element_text(hjust = 0.5),
        axis.title = element_text(size = 17),
        plot.title = element_text(size = 17, face = "bold"))

print(zero_flux_distribution_fold_change)

library(gridExtra)
# Arrange the plots together
combined_plot <- grid.arrange(
    biomass_density, 
    arrangeGrob(zero_flux_distribution, zero_flux_distribution_fold_change, ncol = 1),
    ncol = 2,
    widths = c(2, 1)
)

ggsave("Results/figures/Fig 2.svg", combined_plot, width = 15, height = 8)
