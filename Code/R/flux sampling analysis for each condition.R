setwd('/Users/yunli/GEM-PROSPECT/')
library(tidyverse)
library(ggbeeswarm)
library(FSA)
library(rcompanion)

new_labels <- c("Auto", "Mixo", "Hetero", "Auto_CO2", 
                "Mixo_CO2", "Mixo_hypo", "Mixo_hypo25", "Mixo_hypo75")

# read and process a single condition index
process_condition <- function(i) {
  cond_name <- paste0("Cond", i)
  
  # specific file paths
  f1 <- paste0("Results/flux_sampling/", cond_name, "_sampling.csv")
  f2 <- paste0("Results/flux_sampling/", cond_name, "_sampling_re.csv")
  
  label_name <- new_labels[i]
  
  # read both files and bind them immediately
  bind_rows(read.csv(f1), read.csv(f2)) %>%
    mutate(
      Type = label_name,
      mean_flux = log10(meanFlux) 
    ) %>%
    select(RxnID = RxnIndex, mean_flux, Type)
}

# run the function for 1 to 8 and combine into one dataframe
rxn_mean_flux <- map_dfr(1:8, process_condition)

rxn_mean_flux$Type <- factor(rxn_mean_flux$Type, levels = new_labels)

# Kruskal-Wallis Test and Post-hoc Dunn's Test
# remove non-finite values
clean_data <- rxn_mean_flux %>%
  filter(is.finite(mean_flux))

# Kruskal-Wallis Test
kruskal_res <- kruskal.test(mean_flux ~ Type, data = clean_data)

# run the Post-hoc test (Dunn's Test)
# using Benjamini-Hochberg (bh) adjustment
dunn_res <- dunnTest(mean_flux ~ Type, 
                     data = clean_data, 
                     method = "bh")

# generate the  group letters (a, b, ab, etc.)
cld_res <- cldList(P.adj ~ Comparison,
                   data = dunn_res$res,
                   threshold = 0.05)

# let the letter sit above the maximum value (including outliers) of each group
annotation_df <- clean_data %>%
  group_by(Type) %>%
  summarise(max_y = max(mean_flux, na.rm = TRUE)) %>%
  # join with the letters.
  left_join(cld_res, by = c("Type" = "Group")) %>%
  # add a little buffer to the Y position so the letter isn't touching the dot
  mutate(y_pos = max_y + (max(clean_data$mean_flux) - min(clean_data$mean_flux)) * 0.05)

# plot
ggplot(rxn_mean_flux, aes(x = Type, y = mean_flux)) +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(dodge.width = 0.9, varwidth = TRUE, size = 0.1, alpha = 0.5) +
  geom_text(data = annotation_df, 
            aes(x = Type, y = y_pos, label = Letter), 
            size = 5, 
            vjust = 0,
            color = "black") +
  scale_fill_brewer(palette = "Paired") +
  labs(y = "log10(Mean Flux)") +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = 'transparent', colour = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 15, color = "black"),
    plot.title = element_text(size = 13, face = "bold"),
    axis.title.x = element_blank(),
    legend.position = "right"
)

ggsave("Results/figures/Sup Fig 1.svg", width = 8, height = 5)
