setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)

rxn_list <- read.table("Data/Reactions/list_of_rxns_with_at_least_one_protein_in_both.csv", header = TRUE, sep = ",")

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
max_re <- list(
  omega001 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re_omega_0.01.csv")),
  omega005 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re.csv")),
  omega01 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_Re_omega_0.10.csv"))
)

max_ir <- list(
  omega001 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_omega_0.01.csv")),
  omega005 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8.csv")),
  omega01 = read_and_clean(paste0(topDir, "screens/Max_flux_screen_8_omega_0.10.csv"))
)

sam_re <- list(
  Auto   = read.table(paste0(topDir, "flux_sampling/auto_sampling_re.csv"),   header=TRUE, sep=","),
  Hetero = read.table(paste0(topDir, "flux_sampling/hetero_sampling_re.csv"), header=TRUE, sep=","),
  Mixo   = read.table(paste0(topDir, "flux_sampling/mixo_sampling_re.csv"),   header=TRUE, sep=",")
)

sam_ir <- list(
  Auto   = read.table(paste0(topDir, "flux_sampling/auto_sampling.csv"),     header=TRUE, sep=","),
  Hetero = read.table(paste0(topDir, "flux_sampling/hetero_sampling.csv"),   header=TRUE, sep=","),
  Mixo   = read.table(paste0(topDir, "flux_sampling/mixo_sampling.csv"),     header=TRUE, sep=",")
)

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

results <- lapply(names(max_re), function(al) {
  lapply(2:5, function(th) {
    lapply(names(sam_re), function(tp) {
      df_re <- apply_foldchange_filter(max_re[[al]], sam_re[[tp]], -th)
      df_ir <- apply_foldchange_filter(max_ir[[al]], sam_ir[[tp]], -th)
      cbind(df_re, df_ir[-1])  # drop the EnzymeID column from the irreversibles
    }) %>% setNames(names(sam_re))
  }) %>% setNames(paste0("T", 2:5))
}) %>% setNames(names(max_re))

## recall rates calculation
enzyme_list <- list(
  omega001 = data.frame(EnzymesID = character(), stringsAsFactors = FALSE),
  omega005 = data.frame(EnzymesID = character(),stringsAsFactors = FALSE),
  omega01 = data.frame(EnzymesID = character(),stringsAsFactors = FALSE)
)

# calculate the reaction-gene pairs and find the unique enzyme IDs
ran_gene_pair_count <- list(
  omega001 = 0,
  omega005 = 0,
  omega01 = 0
)

for (i in 1:nrow(rxn_list)) {
  enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
  for (al in names(enzyme_list)) {
    for (j in 1:length(enzymes)) {
      enzyme <- enzymes[j]
      rowindex <- grep(enzyme, results[[al]]$T2$Auto$EnzymeID)
      if (length(rowindex) == 0) next
      enzyme_list[[al]] <- rbind(enzyme_list[[al]], data.frame(EnzymesID = enzyme, stringsAsFactors = FALSE))
      ran_gene_pair_count[[al]]  <- ran_gene_pair_count[[al]] + 1
    } 
  }
}

enzyme_list <- lapply(enzyme_list, unique)

## calculate the recall rate (TP of all the reaction-gene pairs / number of the reaction-gene pairs)
thresholds <- c(2, 3, 4, 5)
types <- c("Auto", "Hetero", "Mixo")
omega <- c("omega001", "omega005", "omega01")

first_case <- rxn_list

# only calculate the gene existence in the gpr rules and mutant library
for (al in omega) {
  for (th in thresholds) {
    for (tp in types) {
      current_data = results[[al]][[paste0("T", th)]][[tp]]
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      for (i in 1:nrow(rxn_list)){
        enzymes <- strsplit(rxn_list$Enzymes[i], ";")[[1]]
        count = 0
        for (j in 1:length(enzymes)) {
          enzyme <- enzymes[j]
          rowindex <- grep(enzyme, current_data$EnzymeID)
          if (length(rowindex) == 0) next
          
          rxn <- rxn_list$RxnIndex[i]
          tmp_count <- sum(!is.na(current_data[[rowindex,rxn]]) & 
                             current_data[[rowindex,rxn]] == 0)
          count <- count + tmp_count
        }
        first_case[[col_name]][i] <- count
      }
    }
  }
}

# calculate the recall rate for the first case
f_recall_rate <- list(
  omega001 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega005 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega01 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE)
)

for (al in omega) {
  for (th in thresholds) {
    for (tp in types) {
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      f_recall_rate[[al]] <- rbind(f_recall_rate[[al]], data.frame(
        Type = tp,
        Threshold = paste0("-", th),
        Recall = sum(first_case[[col_name]]) / ran_gene_pair_count[[al]],
        N_rxns = sum(first_case[[col_name]])))
    }
  }
}

f_recall_rate

## Reaction center point of view
## calculate the recall rate (TP of the group of the reaction-gene pairs by reactions / number of the reactions)
second_case <- first_case

# calculate the recall rate for the second case
s_recall_rate <- list(
  omega001 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega005 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega01 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE)
)

for (al in omega) {
  for (th in thresholds) {
    for (tp in types) {
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      s_recall_rate[[al]] <- rbind(s_recall_rate[[al]], data.frame(
        Type = tp,
        Threshold = paste0("-", th),
        Recall = sum(second_case[[col_name]]>0) / nrow(rxn_list),
        N_rxns = sum(second_case[[col_name]]>0)))
    }
  }
}

s_recall_rate

## Gene center point of view
## calculate the recall rate (TP of the group of the reaction-gene pairs by gene / number of the unique gene)
third_case <- enzyme_list

# only calculate the gene existence in the gpr rules and mutant library
for (al in omega) {
  for (th in thresholds) {
    for (tp in types) {
      current_data = results[[al]][[paste0("T", th)]][[tp]]
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      for (i in 1:nrow(enzyme_list[[al]])){
        count = 0
        
        enzyme <- enzyme_list[[al]]$EnzymesID[i]
        rowindex <- grep(enzyme, current_data$EnzymeID)
        rxnindex <- grep(enzyme, rxn_list$Enzymes)
        
        for (j in 1:length(rxnindex)) {
          rxn <- rxn_list$RxnIndex[rxnindex[j]]
          tmp_count <- sum(!is.na(current_data[[rowindex,rxn]]) &
                             current_data[[rowindex,rxn]] == 0)
          count <- count + tmp_count
        }
        # if the mutant of gene blocks one of the associated reactions, then the count is 1
        if (count > 0) {
          third_case[[al]][[col_name]][i] <- 1
        } else {
          third_case[[al]][[col_name]][i] <- 0
        }
      }
    }
  }
}

# calculate the recall rate for the second case
t_recall_rate <- list(
  omega001 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega005 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE),
  omega01 = data.frame(Type = character(), Threshold = character(), Recall = numeric(), N_rxns = numeric(), stringsAsFactors = FALSE)
)

for (al in omega) {
  for (th in thresholds) {
    for (tp in types) {
      col_name <- paste0(al, "_", tp, "_mutant_T", th)
      t_recall_rate[[al]] <- rbind(t_recall_rate[[al]], data.frame(
        Type = tp,
        Threshold = paste0("-", th),
        Recall = sum(third_case[[al]][[col_name]]) / nrow(enzyme_list[[al]]),
        N_rxns = sum(third_case[[al]][[col_name]])))
    }
  }
}

t_recall_rate

# plot violin plot
library(ggbeeswarm)
library(ggplot2)
library(scales)
plot_recall_vs_count <- function(data, title, y_label_l, y_label_r, x_label, omega_text){
  
  # define the ranges
  primary_range   <- c(0, 0.5)
  secondary_range <- c(0, 70)
  
  rescale_to_secondary <- function(x) {
    scales::rescale(x, from = primary_range, to = secondary_range)
  }
  
  ggplot() +
    # bars: N_rxns
    geom_col(data = data,
             aes(x = Threshold, 
                 y = rescale(N_rxns, from = secondary_range, to = primary_range), 
                 fill = Type),
             position = position_dodge(width = .8),
             width = .7,
             alpha = 0.6) +
    # lines + points
    geom_line(data = data,
              aes(x = Threshold, y = Recall, color = Type, group = Type),
              size = 1) +
    geom_point(data = data,
               aes(x = Threshold, y = Recall, color = Type),
               size = 2) +
    # add omega 
    annotate("text", 
             x = 4.5,
             y = 0.45, 
             label = omega_text, 
             size = 5, 
             color = "black", 
             hjust = 1) +
    # primary and secondary y‐axes
    scale_y_continuous(
      name   = y_label_l,
      limits = primary_range,
      breaks = seq(primary_range[1], primary_range[2], by = 0.1),
      sec.axis = sec_axis(
        ~ rescale_to_secondary(.),
        name   = y_label_r,
        breaks = seq(secondary_range[1], secondary_range[2], by = 10)
      )) +
    scale_fill_manual(name = "Type", values = c(Auto = "#A6CEE3", Hetero = "#1F78B4", Mixo = "#B2DF8A")) +
    scale_color_manual(name = "Type", values = c(Auto = "#A6CEE3", Hetero = "#1F78B4", Mixo = "#B2DF8A")) +
    labs(title = title, x = x_label) +
    theme_bw() +
    theme(
      axis.title.y.left   = element_text(size = 15, color = "black"),
      axis.title.y.right  = element_text(size = 15, color = "black"),
      axis.title.x       = element_text(size = 15, color = "black"),
      axis.text.y.right   = element_text(color = "black"),
      axis.text.y.left    = element_text(color = "black"),
      axis.text         = element_text(size = 12, color = "black"),
      plot.title = element_text(size = 13, face = "bold"),
      legend.position     = "none",
      panel.grid.major.x  = element_blank()
    )
}

# Combine the plots
library(gridExtra)
omega001 <- gridExtra::grid.arrange(
  plot_recall_vs_count(f_recall_rate$omega001, 'Reaction-gene pair perspective', 
                       "Recall Rate", "Number of reaction-gene pairs", "", "omega = 0.01"),
  plot_recall_vs_count(s_recall_rate$omega001, 'Reaction-centric perspective', 
                       "Recall Rate", "Number of reactions", "", "omega = 0.01"),
  plot_recall_vs_count(t_recall_rate$omega001, 'Gene-centric perspective', 
                       "Recall Rate", "Number of genes", "", "omega = 0.01"),
  ncol = 3
)
omega005 <- gridExtra::grid.arrange(
  plot_recall_vs_count(f_recall_rate$omega005, '', 
                       "Recall Rate", "Number of reaction-gene pairs", "", "omega = 0.05"),
  plot_recall_vs_count(s_recall_rate$omega005, '', 
                       "Recall Rate", "Number of reactions", "", "omega = 0.05"),
  plot_recall_vs_count(t_recall_rate$omega005, '', 
                       "Recall Rate", "Number of genes", "", "omega = 0.05"),
  ncol = 3
)
omega01 <- gridExtra::grid.arrange(
  plot_recall_vs_count(f_recall_rate$omega01, '', 
                       "Recall Rate", "Number of reaction-gene pairs", "Threshold", "omega = 0.1"),
  plot_recall_vs_count(s_recall_rate$omega01, '', 
                       "Recall Rate", "Number of reactions", "Threshold", "omega = 0.1"),
  plot_recall_vs_count(t_recall_rate$omega01, '', 
                       "Recall Rate", "Number of genes", "Threshold","omega = 0.1"),
  ncol = 3
)

combined_plot <- gridExtra::grid.arrange(
  omega001,
  omega005,
  omega01,
  nrow = 3)

# Save the combined plot
ggsave("Results/figures/Sup Fig 3.svg", combined_plot, width = 12, height = 13)
