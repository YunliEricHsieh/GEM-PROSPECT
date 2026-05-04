setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(ggplot2)
library(patchwork)

for_validation <- read.table("Results/CataPro/enzymes_and_substrates_for_reference_kcats.csv", header = TRUE, sep = ",")
for_prediction <- read.table("Results/CataPro/enzymes_and_substrates_kcats.csv", header = TRUE, sep = ",")

# Helper function to generate a density plot for a given measure
plot_distribution <- function(val_df, pred_df, pathway, measure_col, xlab_expr, ylab_expr, title_text, annotate_coords) {
    combined <- bind_rows(
    data.frame(Type = "Enzymes in model",   Value = val_df[[measure_col]][val_df$Pathway == pathway]),
    data.frame(Type = "Enzymes candidates", Value = pred_df[[measure_col]][pred_df$Pathway == pathway])
  )
  
  # Calculate means for each group
  means <- combined %>% 
    group_by(Type) %>% 
    summarise(mean_value = mean(Value, na.rm = TRUE), .groups = 'drop')
  
  # Perform the Kolmogorov-Smirnov test
  validate_vals  <- combined$Value[combined$Type == "Enzymes in model"]
  predicted_vals <- combined$Value[combined$Type == "Enzymes candidates"]
  
  ks_res <- ks.test(validate_vals, predicted_vals)
  label  <- paste0("p = ", format(ks_res$p.value, digits = 2))
  
  type_colors <- c("Enzymes in model" = "forestgreen", "Enzymes candidates" = "darkorange")
  
  # Create the density plot
  ggplot(combined, aes(x = Value, y = after_stat(scaled), fill = Type)) +
    geom_density(alpha = 0.5, position = "identity") +
    geom_vline(data = means, aes(xintercept = mean_value, color = Type),
               linetype = "dashed", linewidth = 1, show.legend = FALSE) +
    
    scale_fill_manual(values = type_colors) +
    scale_color_manual(values = type_colors) +
    
    labs(title = title_text, x = xlab_expr, y = ylab_expr, fill = "Type") +
    xlim(-4, 4) +
    
    theme_bw() +
    theme(
      panel.border = element_rect(fill = NA, colour = "black"),
      axis.text    = element_text(size = 12),
      axis.title   = element_text(size = 12),
      legend.text  = element_text(size = 11),
      plot.title   = element_text(size = 10, face = "bold")
    ) +
    annotate("text", x = annotate_coords[1], y = annotate_coords[2], 
             label = label, size = 4, color = "black")
}

# --- Butanoate metabolism ---
TN_butanoate <- plot_distribution(
  for_validation, for_prediction, "Butanoate metabolism", "pred_log10.kcat.s..1..", 
  '', "Density (scaled to 1)", "Butanoate metabolism", c(-2.5, 0.8)
)

MC_butanoate <- plot_distribution(
  for_validation, for_prediction, "Butanoate metabolism", "pred_log10.Km.mM..", 
  '', "Density (scaled to 1)", "", c(2.5, 0.8)
)

# --- Pyrimidine metabolism ---
TN_pyrimidine <- plot_distribution(
  for_validation, for_prediction, "Pyrimidine metabolism", "pred_log10.kcat.s..1..", 
  'log-transformed Kcat', '', "Pyrimidine metabolism", c(-2.5, 0.8)
)

MC_pyrimidine <- plot_distribution(
  for_validation, for_prediction, "Pyrimidine metabolism", "pred_log10.Km.mM..", 
  'log-transformed Km', '', "", c(2.5, 0.8)
)


# --- Glycerolipid metabolism ---
TN_glycerolipid <- plot_distribution(
  for_validation, for_prediction, "Glycerolipid metabolism", "pred_log10.kcat.s..1..", 
  '', '', "Glycerolipid metabolism", c(-2.5, 0.8)
)

MC_glycerolipid <- plot_distribution(
  for_validation, for_prediction, "Glycerolipid metabolism", "pred_log10.Km.mM..", 
  '', '', "", c(2.5, 0.8)
)

final_plot <- guide_area() / 
  (TN_butanoate | TN_pyrimidine | TN_glycerolipid) /
  (MC_butanoate | MC_pyrimidine | MC_glycerolipid) +
  
  plot_layout(
    guides = "collect", 
    heights = c(0.1, 1, 1)
  ) & 
  theme(legend.position = "top")

print(final_plot)

# Save the arranged plot
ggsave("Results/figures/CataPro.svg", plot = final_plot, width = 8, height = 5)
