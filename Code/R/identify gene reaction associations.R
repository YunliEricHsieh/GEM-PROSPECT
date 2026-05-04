setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(tidyr)
library(purrr)

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

# reaction without GPR rules
rxn_list_no_gpr <- read.table("Data/Reactions/list_of_rxns_without_GPR.csv", header = TRUE, sep = ",")

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

# Calculate with Threshold 1 (-1 log10)
results <- lapply(1:1, function(th) {
  lapply(names(sam_re), function(tp) {
    df_re <- apply_foldchange_filter(max_re, sam_re[[tp]], -th)
    df_ir <- apply_foldchange_filter(max_ir, sam_ir[[tp]], -th)
    
    # Safely merge via EnzymeID
    full_join(df_re, df_ir, by = "EnzymeID") 
  }) %>% setNames(names(sam_re))
}) %>% setNames(paste0("T", 1:1))

# Extract datasets for Threshold 1
T1auto   <- results$T1$Auto
T1hetero <- results$T1$Hetero
T1mixo   <- results$T1$Mixo

# Identify Potential Associations (Without GPR)
# Count mutants causing zero flux
rxn_mutants_no_gpr <- rxn_list_no_gpr %>% 
  transmute(RxnIndex, RxnID) %>% 
  mutate(
    Auto_mutant   = map_int(RxnIndex, ~ sum(T1auto[[.x]] == 0, na.rm = TRUE)),
    Hetero_mutant = map_int(RxnIndex, ~ sum(T1hetero[[.x]] == 0, na.rm = TRUE)),
    Mixo_mutant   = map_int(RxnIndex, ~ sum(T1mixo[[.x]] == 0, na.rm = TRUE))
  )

# Filter reactions (0 < mutants < 10) and find intersecting enzymes
potential_association <- rxn_mutants_no_gpr %>%
  filter(
    between(Auto_mutant,   1, 19),
    between(Hetero_mutant, 1, 19),
    between(Mixo_mutant,   1, 19)
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_candidates = {
      rxn <- RxnIndex
      E1 <- T1auto$EnzymeID[!is.na(T1auto[[rxn]]) & T1auto[[rxn]] == 0]
      E2 <- T1hetero$EnzymeID[!is.na(T1hetero[[rxn]]) & T1hetero[[rxn]] == 0]
      E3 <- T1mixo$EnzymeID[!is.na(T1mixo[[rxn]]) & T1mixo[[rxn]] == 0]
      paste(Reduce(intersect, list(E1, E2, E3)), collapse = ";")
    }
  ) %>% ungroup()

# Refine GPR Associations (With GPR)
# Reactions with GPR rules
rxn_list_w_gpr <- read.table("Data/Reactions/list_of_rxns_with_proteins_in_both.csv", header = TRUE, sep = ",")
long_rxn <- rxn_list_w_gpr %>% separate_rows(Enzymes, sep = ";")

blocked <- long_rxn %>% 
  mutate(
    blk_auto   = map2_lgl(Enzymes, RxnIndex, ~{ idx <- match(.x, T1auto$EnzymeID); !is.na(idx) && T1auto[idx, .y] == 0 }),
    blk_hetero = map2_lgl(Enzymes, RxnIndex, ~{ idx <- match(.x, T1hetero$EnzymeID); !is.na(idx) && T1hetero[idx, .y] == 0 }),
    blk_mixo   = map2_lgl(Enzymes, RxnIndex, ~{ idx <- match(.x, T1mixo$EnzymeID); !is.na(idx) && T1mixo[idx, .y] == 0 })
  )

rxn_stats <- blocked %>%
  group_by(RxnIndex, RxnID, Enzymes) %>%
  summarise(
    Auto   = sum(blk_auto, na.rm = TRUE),
    Hetero = sum(blk_hetero, na.rm = TRUE),
    Mixo   = sum(blk_mixo, na.rm = TRUE),
    .groups = 'drop_last'
  ) %>%
  summarise(Auto = sum(Auto), Hetero = sum(Hetero), Mixo = sum(Mixo), .groups = 'drop')

# Calculate Mutants and Intersecting Enzymes
rxn_mutants_w_gpr <- rxn_list_w_gpr %>% 
  transmute(RxnIndex, RxnID) %>% 
  mutate(
    Auto_mutant   = map_int(RxnIndex, ~ sum(T1auto[[.x]] == 0, na.rm = TRUE)),
    Hetero_mutant = map_int(RxnIndex, ~ sum(T1hetero[[.x]] == 0, na.rm = TRUE)),
    Mixo_mutant   = map_int(RxnIndex, ~ sum(T1mixo[[.x]] == 0, na.rm = TRUE))
  )

rxn_list2 <- rxn_list_w_gpr %>%
  left_join(rxn_stats,         by = c("RxnIndex", "RxnID")) %>%
  left_join(rxn_mutants_w_gpr, by = c("RxnIndex", "RxnID"))

potential_association_w_gpr <- rxn_mutants_w_gpr %>%
  filter(
    between(Auto_mutant,   1, 19),
    between(Hetero_mutant, 1, 19),
    between(Mixo_mutant,   1, 19)
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_candidates = {
      rxn <- RxnIndex
      E1 <- T1auto$EnzymeID[!is.na(T1auto[[rxn]]) & T1auto[[rxn]] == 0]
      E2 <- T1hetero$EnzymeID[!is.na(T1hetero[[rxn]]) & T1hetero[[rxn]] == 0]
      E3 <- T1mixo$EnzymeID[!is.na(T1mixo[[rxn]]) & T1mixo[[rxn]] == 0]
      paste(Reduce(intersect, list(E1, E2, E3)), collapse = ";")
    }
  ) %>% ungroup()

# Combine the two result tables
final_association <- potential_association %>%
  select(RxnIndex, RxnID, enzyme_candidates) %>%
  bind_rows(
    potential_association_w_gpr %>% select(RxnIndex, RxnID, enzyme_candidates)
  )

# write the final association table
write.table(
  final_association, 
  file = "Results/potential_gene_reaction_associations.csv", 
  sep = ",", 
  row.names = FALSE, 
  quote = FALSE
)
