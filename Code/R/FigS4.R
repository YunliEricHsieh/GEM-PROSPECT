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

plot_recall_vs_count <- function(data, title, y_label_l, y_label_r, x_label, alpha_text) {
  primary_range   <- c(0, 0.8)
  secondary_range <- c(0, 130)
  rescale_to_secondary <- function(x) scales::rescale(x, from = primary_range, to = secondary_range)
  
  type_colors <- c(Auto = "#A6CEE3", Hetero = "#1F78B4", Mixo = "#B2DF8A")
  
  ggplot(data) +
    # Bars
    geom_col(aes(x = Threshold, y = rescale(N_rxns, from = secondary_range, to = primary_range), fill = Type),
             position = position_dodge(width = 0.8), width = 0.7, alpha = 0.6) +
    # Lines & Points
    geom_line(aes(x = Threshold, y = Recall, color = Type, group = Type), linewidth = 1) +
    geom_point(aes(x = Threshold, y = Recall, color = Type), size = 2) +
    
    # Add Alpha Annotation
    annotate("text", x = 5.5, y = 0.7, label = alpha_text, size = 5, color = "black", hjust = 1) +
    
    # Axes
    scale_y_continuous(
      name = y_label_l, limits = primary_range, breaks = seq(primary_range[1], primary_range[2], by = 0.1),
      sec.axis = sec_axis(~ rescale_to_secondary(.), name = y_label_r, breaks = seq(secondary_range[1], secondary_range[2], by = 10))
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
      axis.text.y.right  = element_text(color = "black"),
      axis.text.y.left   = element_text(color = "black"),
      axis.text          = element_text(size = 12, color = "black"),
      plot.title         = element_text(size = 13, face = "bold"),
      legend.position    = "none",
      panel.grid.major.x = element_blank()
    )
}

rxn_list <- read.table("Data/Reactions/list_of_rxns_with_at_least_one_protein_in_both.csv", header = TRUE, sep = ",")

max_re <- list(
  alpha001 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re_alpha_0.01.csv")),
  alpha005 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re_alpha_0.05.csv")),
  alpha01  = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re.csv"))
)

max_ir <- list(
  alpha001 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_alpha_0.01.csv")),
  alpha005 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_alpha_0.05.csv")),
  alpha01  = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8.csv"))
)

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
results <- lapply(names(max_re), function(al) {
  lapply(1:5, function(th) {
    lapply(names(sam_re), function(tp) {
      df_re <- apply_foldchange_filter(max_re[[al]], sam_re[[tp]], -th)
      df_ir <- apply_foldchange_filter(max_ir[[al]], sam_ir[[tp]], -th)
      full_join(df_re, df_ir, by = "EnzymeID") 
    }) %>% setNames(names(sam_re))
  }) %>% setNames(paste0("T", 1:5))
}) %>% setNames(names(max_re))

#  Recall Rates Calculation
thresholds <- c(1, 2, 3, 4, 5)
types <- c("Auto", "Hetero", "Mixo")
alpha <- c("alpha001", "alpha005", "alpha01")

enzyme_list <- list(
  alpha001 = data.frame(EnzymesID = character(), stringsAsFactors = FALSE),
  alpha005 = data.frame(EnzymesID = character(), stringsAsFactors = FALSE),
  alpha01  = data.frame(EnzymesID = character(), stringsAsFactors = FALSE)
)

ran_gene_pair_count <- list(alpha001 = 0, alpha005 = 0, alpha01 = 0)

for (i in 1:nrow(rxn_list)) {
  enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
  for (al in alpha) {
    for (enzyme in enzymes) {
      rowindex <- grep(enzyme, results[[al]]$T2$Auto$EnzymeID)
      if (length(rowindex) > 0) {
        enzyme_list[[al]] <- rbind(enzyme_list[[al]], data.frame(EnzymesID = enzyme, stringsAsFactors = FALSE))
        ran_gene_pair_count[[al]] <- ran_gene_pair_count[[al]] + 1
      }
    } 
  }
}

enzyme_list <- lapply(enzyme_list, unique)

first_case <- rxn_list

for (al in alpha) {
  for (th in thresholds) {
    for (tp in types) {
      current_data <- results[[al]][[paste0("T", th)]][[tp]]
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      
      for (i in 1:nrow(rxn_list)) {
        enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
        count <- 0
        for (enzyme in enzymes) {
          rowindex <- grep(enzyme, current_data$EnzymeID)
          if (length(rowindex) > 0) {
            rxn <- rxn_list$RxnIndex[i]
            tmp_count <- sum(!is.na(current_data[[rowindex, rxn]]) & current_data[[rowindex, rxn]] == 0)
            count <- count + tmp_count
          }
        }
        first_case[[col_name]][i] <- count
      }
    }
  }
}

# Initialize summary lists
f_recall_rate <- list()
s_recall_rate <- list()
t_recall_rate <- list()

second_case <- first_case
third_case <- enzyme_list

for (al in alpha) {
  # Setup empty dataframes
  f_recall_rate[[al]] <- data.frame(Type=character(), Threshold=character(), Recall=numeric(), N_rxns=numeric(), stringsAsFactors=FALSE)
  s_recall_rate[[al]] <- data.frame(Type=character(), Threshold=character(), Recall=numeric(), N_rxns=numeric(), stringsAsFactors=FALSE)
  t_recall_rate[[al]] <- data.frame(Type=character(), Threshold=character(), Recall=numeric(), N_rxns=numeric(), stringsAsFactors=FALSE)
  
  # Compute Third Case matrix
  for (th in thresholds) {
    for (tp in types) {
      current_data <- results[[al]][[paste0("T", th)]][[tp]]
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      
      for (i in 1:nrow(enzyme_list[[al]])) {
        count <- 0
        enzyme <- enzyme_list[[al]]$EnzymesID[i]
        rowindex <- grep(enzyme, current_data$EnzymeID)
        rxnindex <- grep(enzyme, rxn_list$Enzymes)
        
        for (idx in rxnindex) {
          rxn <- rxn_list$RxnIndex[idx]
          tmp_count <- sum(!is.na(current_data[[rowindex, rxn]]) & current_data[[rowindex, rxn]] == 0)
          count <- count + tmp_count
        }
        third_case[[al]][[col_name]][i] <- ifelse(count > 0, 1, 0)
      }
      
      # Extract summaries for all 3 cases
      thresh_label <- paste0("-", th)
      
      f_recall_rate[[al]] <- rbind(f_recall_rate[[al]], data.frame(
        Type = tp, Threshold = thresh_label, 
        Recall = sum(first_case[[col_name]]) / ran_gene_pair_count[[al]], 
        N_rxns = sum(first_case[[col_name]])
      ))
      
      s_recall_rate[[al]] <- rbind(s_recall_rate[[al]], data.frame(
        Type = tp, Threshold = thresh_label, 
        Recall = sum(second_case[[col_name]] > 0) / nrow(rxn_list), 
        N_rxns = sum(second_case[[col_name]] > 0)
      ))
      
      t_recall_rate[[al]] <- rbind(t_recall_rate[[al]], data.frame(
        Type = tp, Threshold = thresh_label, 
        Recall = sum(third_case[[al]][[col_name]]) / nrow(enzyme_list[[al]]), 
        N_rxns = sum(third_case[[al]][[col_name]])
      ))
    }
  }
}

# Generate plots for alpha 0.01
a1 <- plot_recall_vs_count(f_recall_rate$alpha001, 'Reaction-gene pair perspective', "Recall Rate", "Number of reaction-gene pairs", "", "alpha = 0.01")
a2 <- plot_recall_vs_count(s_recall_rate$alpha001, 'Reaction-centric perspective',   "Recall Rate", "Number of reactions",           "", "alpha = 0.01")
a3 <- plot_recall_vs_count(t_recall_rate$alpha001, 'Gene-centric perspective',       "Recall Rate", "Number of genes",               "", "alpha = 0.01")

# Generate plots for alpha 0.05
b1 <- plot_recall_vs_count(f_recall_rate$alpha005, '', "Recall Rate", "Number of reaction-gene pairs", "", "alpha = 0.05")
b2 <- plot_recall_vs_count(s_recall_rate$alpha005, '', "Recall Rate", "Number of reactions",           "", "alpha = 0.05")
b3 <- plot_recall_vs_count(t_recall_rate$alpha005, '', "Recall Rate", "Number of genes",               "", "alpha = 0.05")

# Generate plots for alpha 0.1
c1 <- plot_recall_vs_count(f_recall_rate$alpha01, '', "Recall Rate", "Number of reaction-gene pairs", "Threshold", "alpha = 0.1")
c2 <- plot_recall_vs_count(s_recall_rate$alpha01, '', "Recall Rate", "Number of reactions",           "Threshold", "alpha = 0.1")
c3 <- plot_recall_vs_count(t_recall_rate$alpha01, '', "Recall Rate", "Number of genes",               "Threshold", "alpha = 0.1")

# Combine all 9 plots into a 3x3 grid using Patchwork
combined_plot <- (a1 | a2 | a3) / 
  (b1 | b2 | b3) / 
  (c1 | c2 | c3)

# Save the combined plot
ggsave("Results/figures/Sup Fig 4.svg", combined_plot, width = 12, height = 13)
