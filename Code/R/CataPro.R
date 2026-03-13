setwd('/Users/yunli/GEM-PROSPECT/')
library(ggplot2)
library(dplyr)

for_validation <- read.table("Results/CataPro/enzymes_and_substrates_for_reference_kcats.csv", header = TRUE, sep = ",")
for_prediction <- read.table("Results/CataPro/enzymes_and_substrates_kcats.csv", header = TRUE, sep = ",")

# Helper function to generate a density plot for a given measure
plot_distribution <- function(val_df, pred_df, pathway, measure_col, xlab_expr, ylab_expr, title_text, annotate_coords) {
  # Combine data for the given measure from validation and prediction files
  combined <- rbind(
    data.frame(Type = "Enzymes in model", Value = val_df[[measure_col]][val_df$Pathway == pathway]),
    data.frame(Type = "Enzymes candidates", Value = pred_df[[measure_col]][pred_df$Pathway == pathway])
  )
  
  # Calculate means for each group
  means <- combined %>% 
    group_by(Type) %>% 
    summarise(mean_value = mean(Value, na.rm = TRUE)) %>% 
    ungroup()
  
  # Perform the Kolmogorov-Smirnov test
  validate_vals <- combined$Value[combined$Type == "Enzymes in model"]
  predicted_vals <- combined$Value[combined$Type == "Enzymes candidates"]
  ks_res <- ks.test(validate_vals, predicted_vals)
  label <- paste0("p = ", format(ks_res$p.value, digits = 2))
  
  # Create the density plot with vertical mean lines and annotation
  ggplot(combined, aes(x = Value,after_stat(scaled), fill = Type)) +
    geom_density(alpha = 0.5, position = "identity") +
    geom_vline(data = means, aes(xintercept = mean_value, color = Type),
               linetype = "dashed", size = 1, show.legend = FALSE) +
    scale_fill_manual(values = c("Enzymes in model" = "forestgreen", "Enzymes candidates" = "darkorange")) +
    scale_color_manual(values = c("Enzymes in model" = "forestgreen", "Enzymes candidates" = "darkorange")) +
    labs(title = title_text,
         x = xlab_expr,
         y = ylab_expr,
         fill = "Type") +
    xlim(-3, 3) +
    theme_bw() +
    theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 11),
          plot.title = element_text(size = 10, face = "bold")) +
    annotate("text", x = annotate_coords[1], y = annotate_coords[2], 
             label = label, size = 3, color = "black")
}

# Generate plots for each measure
# Purine metabolism
# Turnover (kcat)
TN_purine <- plot_distribution(for_validation, for_prediction, "Purine metabolism",
                          "pred_log10.kcat.s..1..", 
                          '',
                          "Density (scaled to 1)",
                          "Purine metabolism", 
                          c(-1.5, 0.75))
print(TN_purine)

# Michaelis constant (Km)
MC_purine <- plot_distribution(for_validation, for_prediction, "Purine metabolism",
                          "pred_log10.Km.mM..", 
                          '',
                          "Density (scaled to 1)",
                          "", 
                          c(1.7, 0.75))
print(MC_purine)


# Pyruvate metabolism
# Turnover (kcat)
TN_pyruvate <- plot_distribution(for_validation, for_prediction, "Pyruvate metabolism",
                          "pred_log10.kcat.s..1..", 
                          'log-transformed Kcat',
                          '',
                          "Pyruvate metabolism", 
                          c(-1, 0.75))
print(TN_pyruvate)

# Michaelis constant (Km)
MC_pyruvate <- plot_distribution(for_validation, for_prediction, "Pyruvate metabolism",
                          "pred_log10.Km.mM..", 
                          'log-transformed Km',
                          '',
                          "", 
                          c(1.7, 0.8))
print(MC_pyruvate)

# Urea degradation
# Turnover (kcat)
TN_urea <- plot_distribution(for_validation, for_prediction, "Urea degradation",
                          "pred_log10.kcat.s..1..", 
                          '',
                          '',
                          "Urea degradation", 
                          c(-1.5, 0.85))
print(TN_urea)

# Michaelis constant (Km)
MC_urea <- plot_distribution(for_validation, for_prediction, "Urea degradation",
                          "pred_log10.Km.mM..", 
                          '',
                          '',
                          "", 
                          c(2, 0.85))
print(MC_urea)

# Glycerolipid metabolism
# Turnover (kcat)
reference <- data.frame(Type = "Enzymes in model", 
             Value = for_validation[["pred_log10.kcat.s..1.."]][for_validation$Pathway == "Glycerolipid metabolism"])

prediction <- data.frame(Type = "Enzymes in model", 
                        Value = for_prediction[["pred_log10.kcat.s..1.."]][for_prediction$Pathway == "Glycerolipid metabolism"])

# Calculate means for each group
means <- reference %>% 
  group_by(Type) %>% 
  summarise(mean_value = mean(Value, na.rm = TRUE)) %>% 
  ungroup()

# Create the density plot with vertical mean lines and annotation
TN_glycerolipid <- ggplot(reference, aes(x = Value, after_stat(scaled), fill = Type)) +
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(data = means, aes(xintercept = mean_value, color = Type),
             linetype = "dashed", size = 1, show.legend = FALSE) +
  geom_point(
    data  = prediction,
    aes(x = Value, y = 0),
    shape = 17,
    size  = 3,
    color = "darkorange",
    position = position_jitter(height = 0.05)) +
  scale_fill_manual(values = c("Enzymes in model" = "forestgreen")) +
  scale_color_manual(values = c("Enzymes in model" = "forestgreen")) +
  labs(title = "Glycerolipid metabolism",
       x = '',
       y = "",
       fill = "Type") +
  xlim(-3, 3) +
  theme_bw() +
  theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 10, face = "bold")) 

print(TN_glycerolipid)

# Michaelis constant (Km)
reference <- data.frame(Type = "Enzymes in model", 
                        Value = for_validation[["pred_log10.Km.mM.."]][for_validation$Pathway == "Glycerolipid metabolism"])

prediction <- data.frame(Type = "Enzymes in model", 
                         Value = for_prediction[["pred_log10.Km.mM.."]][for_prediction$Pathway == "Glycerolipid metabolism"])

# Calculate means for each group
means <- reference %>% 
  group_by(Type) %>% 
  summarise(mean_value = mean(Value, na.rm = TRUE)) %>% 
  ungroup()

# Create the density plot with vertical mean lines and annotation
MC_glycerolipid <- ggplot(reference, aes(x = Value, after_stat(scaled), fill = Type)) +
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(data = means, aes(xintercept = mean_value, color = Type),
             linetype = "dashed", size = 1, show.legend = FALSE) +
  geom_point(
    data  = prediction,
    aes(x = Value, y = 0),
    shape = 17,
    size  = 3,
    color = "darkorange",
    position = position_jitter(height = 0.05)) +
  scale_fill_manual(values = c("Enzymes in model" = "forestgreen")) +
  scale_color_manual(values = c("Enzymes in model" = "forestgreen")) +
  labs(title = "",
       x = '',
       y = "",
       fill = "Type") +
  xlim(-3, 3) +
  theme_bw() +
  theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        plot.title = element_text(size = 10, face = "bold")) 

print(MC_glycerolipid)

# Arrange the figures in the specified order and layout
library(ggpubr)
final_plot <- ggarrange(
  ggarrange(TN_purine, MC_purine, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top"),
  ggarrange(TN_pyruvate, MC_pyruvate, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top"),
  ggarrange(TN_urea, MC_urea, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top"),
  ggarrange(TN_glycerolipid, MC_glycerolipid, ncol = 1, nrow = 2, common.legend = TRUE, legend = "top"),
  ncol = 4, nrow = 1, common.legend = TRUE, legend = "top"
)

print(final_plot)

# Save the arranged plot
ggsave("Results/figures/Fig 5.svg", plot = final_plot, width = 8, height = 5)
