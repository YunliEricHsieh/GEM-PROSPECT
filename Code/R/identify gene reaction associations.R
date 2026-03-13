setwd('/Users/yunli/GEM-PROSPECT/')

# identify potential gene-reaction associations
# reaction without GPR rules
rxn_list <- read.table("Data/Reactions/list_of_rxns_without_GPR.csv", header = TRUE, sep = ",")

# read & clean one table function
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

# log10 fold-change analysis function
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

# fold change analysis
library(dplyr)
results <- lapply(2:2, function(th) {
  lapply(names(sam_re), function(tp) {
    df_re <- apply_foldchange_filter(max_re, sam_re[[tp]], -th)
    df_ir <- apply_foldchange_filter(max_ir, sam_ir[[tp]], -th)
    cbind(df_re, df_ir[-1])  # drop the EnzymeID column from the irreversibles
  }) %>% setNames(names(sam_re))
}) %>% setNames(paste0("T", 2:2))

T2auto   <- results$T2$Auto
T2hetero <- results$T2$Hetero
T2mixo   <- results$T2$Mixo

# identify reactions with zero flux
library(purrr)
rxn_list <- rxn_list %>% 
  transmute(RxnIndex, RxnID) %>% 
  mutate(
    Auto_mutant   = map_int(RxnIndex, ~ sum(T2auto[[.x]] == 0, na.rm = TRUE)),
    Hetero_mutant = map_int(RxnIndex, ~ sum(T2hetero[[.x]] == 0, na.rm = TRUE)),
    Mixo_mutant   = map_int(RxnIndex, ~ sum(T2mixo[[.x]] == 0, na.rm = TRUE))
  )

# find the rxns have reasonable associated enzyme (n < 20)
potential_association <- rxn_list %>%
  # keep only the reactions with 0 < mutants < 20 in every condition
  filter(
    between(Auto_mutant,   1, 19),
    between(Hetero_mutant, 1, 19),
    between(Mixo_mutant,   1, 19)
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_candidates = {
      rxn <- RxnIndex
      E1   <- T2auto$EnzymeID[!is.na(T2auto[[rxn]]) & T2auto[[rxn]] == 0]
      E2 <- T2hetero$EnzymeID[!is.na(T2hetero[[rxn]]) & T2hetero[[rxn]] == 0]
      E3  <- T2mixo$EnzymeID[!is.na(T2mixo[[rxn]]) & T2mixo[[rxn]] == 0]
      paste(Reduce(intersect, list(E1, E2, E3)), collapse = ";")
    }
  )

## refine GPR associations
# reaction with GPR rules
rxn_list1 <- read.table("Data/Reactions/list_of_rxns_with_proteins_in_both.csv", header = TRUE, sep = ",")

library(tidyr)
long_rxn <- rxn_list1 %>% 
  separate_rows(Enzymes, sep = ";")

# find the reactions with zero flux 
blocked <- long_rxn %>% 
  mutate(
    blk_auto   = map2_lgl(Enzymes, RxnIndex, ~{
      idx <- match(.x, T2auto$EnzymeID)
      !is.na(idx) && T2auto[idx, .y] == 0
    }),
    blk_hetero = map2_lgl(Enzymes, RxnIndex, ~{
      idx <- match(.x, T2hetero$EnzymeID)
      !is.na(idx) && T2hetero[idx, .y] == 0
    }),
    blk_mixo   = map2_lgl(Enzymes, RxnIndex, ~{
      idx <- match(.x, T2mixo$EnzymeID)
      !is.na(idx) && T2mixo[idx, .y] == 0
    })
  )

rxn_stats <- blocked %>%
  group_by(RxnIndex, RxnID, Enzymes) %>%
  summarise(
    Auto   = sum(blk_auto, na.rm = TRUE),
    Hetero = sum(blk_hetero, na.rm = TRUE),
    Mixo   = sum(blk_mixo, na.rm = TRUE),
    .groups = 'drop_last'
  ) %>%
  summarise(
    Auto   = sum(Auto),
    Hetero = sum(Hetero),
    Mixo   = sum(Mixo),
    .groups = 'drop'
  )

rxn_mutants <- rxn_list1 %>% 
  transmute(RxnIndex, RxnID) %>% 
  mutate(
    Auto_mutant   = map_int(RxnIndex, ~ sum(T2auto[[.x]] == 0, na.rm = TRUE)),
    Hetero_mutant = map_int(RxnIndex, ~ sum(T2hetero[[.x]] == 0, na.rm = TRUE)),
    Mixo_mutant   = map_int(RxnIndex, ~ sum(T2mixo[[.x]] == 0, na.rm = TRUE))
  )

rxn_list2 <- rxn_list1 %>%
  left_join(rxn_stats,   by = c("RxnIndex", "RxnID")) %>%
  left_join(rxn_mutants, by = c("RxnIndex", "RxnID"))

potential_association_w_gpr <- rxn_mutants %>%
  filter(
    between(Auto_mutant,   1, 19),
    between(Hetero_mutant, 1, 19),
    between(Mixo_mutant,   1, 19),
  ) %>%
  rowwise() %>%
  mutate(
    enzyme_candidates = {
      rxn <- RxnIndex
      E1   <- T2auto$EnzymeID[!is.na(T2auto[[rxn]]) & T2auto[[rxn]] == 0]
      E2 <- T2hetero$EnzymeID[!is.na(T2hetero[[rxn]]) & T2hetero[[rxn]] == 0]
      E3  <- T2mixo$EnzymeID[!is.na(T2mixo[[rxn]]) & T2mixo[[rxn]] == 0]
      paste(Reduce(intersect, list(E1, E2, E3)), collapse = ";")
    }
  )

# combine two result tables
final_association <- potential_association %>%
  select(RxnIndex, RxnID, enzyme_candidates) %>%
  bind_rows(potential_association_w_gpr %>%
              select(RxnIndex, RxnID, enzyme_candidates))

# write the final association table
write.table(final_association, 
            file = "Results/potential_gene_reaction_associations.csv", 
            sep = ",", row.names = FALSE, quote = FALSE)
