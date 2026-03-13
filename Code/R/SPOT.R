# SPOT prediction score distribution
setwd('/Users/yunli/GEM-PROSPECT/')

library(dplyr)

#  read excel file
library(readxl)
validate_data <- read_excel("Results/SPOT/SPOT_prediction_validate.xlsx")
predicted_data <- read_excel("Results/SPOT/candidate protein of transporter.xlsx")

#  combine the prediction score for distribution plot
combined_scores <- rbind(
  data.frame(Type = "Transport proteins in model", Score = validate_data$`Prediction score`),
  data.frame(Type = "Candidate of transport proteins", Score = predicted_data$`Prediction score`)
)

# Compute means for each Source
means <- combined_scores %>% 
  group_by(Type) %>% 
  summarise(mean_score = mean(Score, na.rm = TRUE))

# find the mode value
mode_value <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode_value(validate_data$`Prediction score`)
mode_value(predicted_data$`Prediction score`)

# Plot distribution and add vertical lines at the mean values
SPOT<-ggplot(combined_scores, aes(x = Score, fill = Type)) +
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(data = means, aes(xintercept = mean_score, color = Type),
             linetype = "dashed", size = 1, show.legend = FALSE) +
  scale_fill_manual(values = c("Transport proteins in model" = "forestgreen", "Candidate of transport proteins" = "darkorange")) +
  scale_color_manual(values = c("Transport proteins in model" = "forestgreen", "Candidate of transport proteins" = "darkorange")) +
  labs(title = "",
       x = "Prediction Score",
       y = "Density",
       fill = "")+
  theme_bw() +
  theme(panel.border = element_rect(fill = 'transparent', colour = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1))

print(SPOT)