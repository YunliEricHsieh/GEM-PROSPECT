setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(stringr)

measure <- read.csv("Data/Mutant_phenotypes_table_filtered_with_uniportIDs.csv",
                    stringsAsFactors = FALSE) %>%
  select(-Feature, -Mutant_ID, -Confidence_level, -Gene_Name)

auto.hetero = rowMeans(measure[, c("c_R5_S_233_TP_Air_I_v_R5_S_217_Control_TAP_Dark_II",
                                   "c_R5_S_234_TP_Air_II_v_R5_S_217_Control_TAP_Dark_II",
                                   "c_R5_S_235_TP_Air_III_v_R5_S_217_Control_TAP_Dark_II")], na.rm = TRUE)

auto.mixo = rowMeans(measure[, c("c_R3_S_008_TP_Air_II_v_R3_A_002_Control_TAP_Light_average",
                                  "c_R3_S_007_TP_Air_I_v_R3_A_002_Control_TAP_Light_average")], na.rm = TRUE)

auto.auto_CO2 = rowMeans(measure[, c("c_R5_S_233_TP_Air_I_v_R5_A_012_TP_CO2_average",
                                     "c_R3_S_008_TP_Air_II_v_R3_A_011_TP_CO2_average",
                                     "c_R5_S_234_TP_Air_II_v_R5_A_012_TP_CO2_average",
                                     "c_R3_S_072_TP_Air_v_R3_S_071_TP_CO2",
                                     "c_R3_S_007_TP_Air_I_v_R3_A_011_TP_CO2_average",
                                  "c_R5_S_235_TP_Air_III_v_R5_A_012_TP_CO2_average")], na.rm = TRUE)

mixo.mixo_CO2 = rowMeans(measure[, c("c_R3_S_001_Control_TAP_Light_I_v_R3_S_011_TAP_Light_CO2",
                                  "c_R3_S_002_Control_TAP_Light_II_v_R3_S_011_TAP_Light_CO2")], na.rm = TRUE)

mixo.hetero = rowMeans(measure[, c("c_R3_S_003_Control_TAP_Light_III_v_R3_A_006_Control_TAP_Dark_av",
                                   "c_R3_S_002_Control_TAP_Light_II_v_R3_A_006_Control_TAP_Dark_ave",
                                   "c_R3_S_063_Control_TAP_Light_v_R3_A_006_Control_TAP_Dark_averag",
                                     "c_R3_S_001_Control_TAP_Light_I_v_R3_A_006_Control_TAP_Dark_aver")], na.rm = TRUE)

mixo_NaCl.mixo = measure$c_R3_S_050_NaCl_50_v_R3_A_002_Control_TAP_Light_average

mixo_P.mixo = measure$c_R3_S_033_Low_P_0_01_mM_v_R3_A_002_Control_TAP_Light_average

mixo_N.mixo = measure$c_R3_S_032_Low_N_0_7_mM_v_R3_A_002_Control_TAP_Light_average

mixo_hypo10.mixo = measure$c_R3_S_060_Hypoosmotic_10_v_R3_A_002_Control_TAP_Light_average

mixo_hypo25.mixo = measure$c_R4_S_146_Hypoosmotic_25_v_R4_A_003_Control_TAP_Light_average

mixo_hypo75.mixo = rowMeans(measure[, c("c_R4_S_147_Hypoosmotic_75_v_R4_A_003_Control_TAP_Light_average",
                                   "c_R3_S_062_Hypoosmotic_75_v_R3_A_002_Control_TAP_Light_average")], na.rm = TRUE)

GeneID <- unique(measure$Gene)

UniProtID <- numeric()
auto_hetero <- numeric()
auto_mixo <- numeric()
auto_auto_CO2 <- numeric()
mixo_mixo_CO2 <- numeric()
mixo_hetero <- numeric()
mixo_NaCl_mixo <- numeric()
mixo_P_mixo <- numeric()
mixo_N_mixo <- numeric()
mixo_hypo10_mixo <- numeric()
mixo_hypo25_mixo <- numeric()
mixo_hypo75_mixo <- numeric()

for (i in 1:length(GeneID)){
  UniProtID[i] <- measure$UniProtID[which(measure$Gene == GeneID[i])]
  auto_hetero[i] <- mean(auto.hetero[which(measure$Gene == GeneID[i])],na.rm = T)
  auto_mixo[i] <- mean(auto.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  auto_auto_CO2[i] <- mean(auto.auto_CO2[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_mixo_CO2[i] <- mean(mixo.mixo_CO2[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_hetero[i] <- mean(mixo.hetero[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_NaCl_mixo[i] <- mean(mixo_NaCl.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_P_mixo[i] <- mean(mixo_P.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_N_mixo[i] <- mean(mixo_N.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_hypo10_mixo[i] <- mean(mixo_hypo10.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_hypo25_mixo[i] <- mean(mixo_hypo25.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
  mixo_hypo75_mixo[i] <- mean(mixo_hypo75.mixo[which(measure$Gene == GeneID[i])],na.rm = T)
}

new_table <- as.data.frame(cbind(GeneID, UniProtID, auto_hetero, auto_mixo, auto_auto_CO2,
                                 mixo_mixo_CO2, mixo_hetero, mixo_NaCl_mixo, mixo_P_mixo,
                                 mixo_N_mixo, mixo_hypo10_mixo, mixo_hypo25_mixo, mixo_hypo75_mixo))

new_table <- new_table[!grepl("&", new_table$GeneID), ]

write.csv(new_table, 'Data/Mutant_phenotypes_table_filtered_final.csv', row.names = F)
