setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
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

plot_recall_vs_count <- function(data, title, y_label_l, y_label_r, x_label){
  primary_range   <- c(0, 0.5)
  secondary_range <- c(0, 10)
  
  rescale_to_secondary <- function(x) scales::rescale(x, from = primary_range, to = secondary_range)
  
  # Define colors once
  type_colors <- c(Auto = "#A6CEE3", Hetero = "#1F78B4", Mixo = "#B2DF8A")
  
  ggplot(data) +
    # Bars (Secondary Axis)
    geom_col(aes(x = Threshold, y = rescale(N_rxns, from = secondary_range, to = primary_range), fill = Type),
             position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
    # Lines & Points (Primary Axis)
    geom_line(aes(x = Threshold, y = Recall, color = Type, group = Type), size = 1) +
    geom_point(aes(x = Threshold, y = Recall, color = Type), size = 2) +
    
    # Axes and Scales
    scale_y_continuous(
      name = y_label_l, limits = primary_range, breaks = seq(primary_range[1], primary_range[2], by = 0.1),
      sec.axis = sec_axis(~ rescale_to_secondary(.), name = y_label_r, breaks = seq(secondary_range[1], secondary_range[2], by = 1))
    ) +
    scale_fill_manual(name = "Type", values = type_colors) +
    scale_color_manual(name = "Type", values = type_colors) +
    labs(title = title, x = x_label) +
    
    # Theme
    theme_bw() +
    theme(
      axis.title.y.left  = element_text(size = 15, color = "black"),
      axis.title.y.right = element_text(size = 15, color = "black"),
      axis.title.x       = element_text(size = 15, color = "black"),
      axis.text          = element_text(size = 12, color = "black"),
      plot.title         = element_text(size = 13, face = "bold"),
      legend.position    = "none",
      panel.grid.major.x = element_blank()
    )
}

rxn_list <- read.table("Data/Reactions/list_of_rxns_with_at_least_one_protein_in_both.csv", header = TRUE, sep = ",")
# Remove reactions associated with multiple enzymes (ones containing a semicolon)
rxn_list <- rxn_list[!grepl(";", rxn_list$Enzymes), ]

max_re <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re.csv"))
max_ir <- read_and_clean(paste0(topDir, "screens/Max_flux_screen_8.csv"))

sam_re <- list(
  Auto   = read.table(paste0(topDir, "flux_sampling/auto_sampling_re.csv"),   header = TRUE, sep = ","),
  Hetero = read.table(paste0(topDir, "flux_sampling/hetero_sampling_re.csv"), header = TRUE, sep = ","),
  Mixo   = read.table(paste0(topDir, "flux_sampling/mixo_sampling_re.csv"),   header = TRUE, sep = ",")
)

sam_ir <- list(
  Auto   = read.table(paste0(topDir, "flux_sampling/auto_sampling.csv"),      header = TRUE, sep = ","),
  Hetero = read.table(paste0(topDir, "flux_sampling/hetero_sampling.csv"),    header = TRUE, sep = ","),
  Mixo   = read.table(paste0(topDir, "flux_sampling/mixo_sampling.csv"),      header = TRUE, sep = ",")
)

# Fold-change analysis
results <- lapply(1:5, function(th) {
  lapply(names(sam_re), function(tp) {
    df_re <- apply_foldchange_filter(max_re, sam_re[[tp]], -th)
    df_ir <- apply_foldchange_filter(max_ir, sam_ir[[tp]], -th)
    full_join(df_re, df_ir, by = "EnzymeID") 
  }) %>% setNames(names(sam_re))
}) %>% setNames(paste0("T", 1:5))

#  Recall Rates Calculation
thresholds <- c(1, 2, 3, 4, 5)
types <- c("Auto", "Hetero", "Mixo")

# Calculate Valid Enzyme List
enzyme_list <- data.frame(EnzymesID = character(), stringsAsFactors = FALSE)
ran_gene_pair_count <- 0 

for (i in 1:nrow(rxn_list)) {
  enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
  for (enzyme in enzymes) {
    if (length(grep(enzyme, results$T2$Auto$EnzymeID)) > 0) {
      enzyme_list <- rbind(enzyme_list, data.frame(EnzymesID = enzyme, stringsAsFactors = FALSE))
      ran_gene_pair_count <- ran_gene_pair_count + 1
    }
  } 
}
enzyme_list <- unique(enzyme_list)

# --- Reaction-Gene Pair Perspective (First Case) ---
first_case <- rxn_list

for (th in thresholds) {
  for (tp in types) {
    current_data <- results[[paste0("T", th)]][[tp]]
    col_name <- paste0(tp, "_mutant_T", th)
    
    for (i in 1:nrow(rxn_list)) {
      enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
      count <- 0
      
      for (enzyme in enzymes) {
        rowindex <- grep(enzyme, current_data$EnzymeID)
        if (length(rowindex) > 0) {
          rxn <- rxn_list$RxnIndex[i]
          # Suppress NA errors during count
          tmp_count <- sum(!is.na(current_data[[rowindex, rxn]]) & current_data[[rowindex, rxn]] == 0)
          count <- count + tmp_count
        }
      }
      first_case[[col_name]][i] <- count
    }
  }
}

# Extract Summary
f_recall_rate <- expand.grid(Type = types, Threshold = thresholds, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    col_name = paste0(Type, "_mutant_T", Threshold),
    Recall = sum(first_case[[col_name]]) / ran_gene_pair_count,
    N_rxns = sum(first_case[[col_name]]),
    Threshold = paste0("-", Threshold)
  ) %>% ungroup() %>% select(-col_name)

# --- Reaction-Centric Perspective (Second Case) ---
second_case <- first_case

s_recall_rate <- expand.grid(Type = types, Threshold = thresholds, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    col_name = paste0(Type, "_mutant_T", Threshold),
    Recall = sum(second_case[[col_name]] > 0) / nrow(rxn_list),
    N_rxns = sum(second_case[[col_name]] > 0),
    Threshold = paste0("-", Threshold)
  ) %>% ungroup() %>% select(-col_name)

# --- Gene-Centric Perspective (Third Case) ---
third_case <- enzyme_list

for (th in thresholds) {
  for (tp in types) {
    current_data <- results[[paste0("T", th)]][[tp]]
    col_name <- paste0(tp, "_mutant_T", th)
    
    for (i in 1:nrow(enzyme_list)) {
      enzyme <- enzyme_list$EnzymesID[i]
      rowindex <- grep(enzyme, current_data$EnzymeID)
      rxnindex <- grep(enzyme, rxn_list$Enzymes)
      
      count <- 0
      for (idx in rxnindex) {
        rxn <- rxn_list$RxnIndex[idx]
        tmp_count <- sum(!is.na(current_data[[rowindex, rxn]]) & current_data[[rowindex, rxn]] == 0)
        count <- count + tmp_count
      }
      third_case[[col_name]][i] <- ifelse(count > 0, 1, 0)
    }
  }
}

t_recall_rate <- expand.grid(Type = types, Threshold = thresholds, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    col_name = paste0(Type, "_mutant_T", Threshold),
    Recall = sum(third_case[[col_name]]) / nrow(enzyme_list),
    N_rxns = sum(third_case[[col_name]]),
    Threshold = paste0("-", Threshold)
  ) %>% ungroup() %>% select(-col_name)


p_f <- plot_recall_vs_count(f_recall_rate, 'Reaction-gene pair perspective', 
                            "Recall Rate", "Number of reaction-gene pairs", "Threshold")
p_s <- plot_recall_vs_count(s_recall_rate, 'Reaction-centric perspective', 
                            "Recall Rate", "Number of reactions", "Threshold")
p_t <- plot_recall_vs_count(t_recall_rate, 'Gene-centric perspective', 
                            "Recall Rate", "Number of genes", "Threshold")

combined_plot <- p_f | p_s | p_t

# Save the combined plot
ggsave("Results/figures/Sup Fig 2.svg", combined_plot, width = 12, height = 4)
