# SPOT prediction score distribution
setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)
library(ggplot2)
library(readxl)

#  read excel file
validate_data  <- read_excel("Results/SPOT/SPOT_prediction_validate.xlsx")
predicted_data <- read_excel("Results/SPOT/potential_transporter_SPOT.xlsx")

#  combine the prediction score for distribution plot
combined_scores <- bind_rows(
  data.frame(Type = "Transport proteins in model",   Score = validate_data$`Prediction score`),
  data.frame(Type = "Candidate of transport proteins", Score = predicted_data$`Prediction score`)
)

# Compute means for each Source
means <- combined_scores %>% 
  group_by(Type) %>% 
  summarise(mean_score = mean(Score, na.rm = TRUE), .groups = 'drop')

# find the mode value
get_mode <- function(x) {
  x <- na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode_validate <- get_mode(validate_data$`Prediction score`)
mode_predicted <- get_mode(predicted_data$`Prediction score`)

cat("Mode (Transport proteins in model):", mode_validate, "\n")
cat("Mode (Candidate of transport proteins):", mode_predicted, "\n")

# Plot distribution and add vertical lines at the mean values
type_colors <- c(
  "Transport proteins in model"     = "forestgreen", 
  "Candidate of transport proteins" = "darkorange"
)

plot_SPOT <- ggplot(combined_scores, aes(x = Score, fill = Type)) +
  geom_density(alpha = 0.5, position = "identity") +
  
  # Note: `size` is deprecated in newer ggplot2 versions for lines; `linewidth` is correct
  geom_vline(data = means, aes(xintercept = mean_score, color = Type),
             linetype = "dashed", linewidth = 1, show.legend = FALSE) +
  
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  
  labs(
    x = "Prediction Score",
    y = "Density",
    fill = NULL # Removes the legend title directly
  ) +
  
  theme_bw() +
  theme(
    # Changed fill = 'transparent' to fill = NA (standard ggplot syntax for hollow rects)
    panel.border = element_rect(fill = NA, colour = "black"), 
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    
    # Legend formatting inside the plot
    legend.text = element_text(size = 12),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    # Optional: Adds a subtle white background to the legend so plot lines don't run through the text
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3) 
  )

print(plot_SPOT)

# Save the plot
ggsave("Results/figures/SPOT.svg", plot = final_plot, width = 8, height = 5)