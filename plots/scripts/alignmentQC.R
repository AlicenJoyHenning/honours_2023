# Alignment Quality 

##### Loading dependencies ####
library(ggplot2)
install.packages("extrafont")
library(extrafont)
font_import(pattern = "EBGaramond")

# BLUE "#a0d9e9"

##### Merged Plot (all datasets) ####
Datasets <- c(rep("alpha" , 2), rep("lambda", 2), rep("untreated", 2)) 
Metrics <- rep(c("Alignment Rate" , "Percentage Whitelist" ),3)
Percentage <- c(56.6, 95.6, 55.3, 95.7, 54.1, 96.1)

alignQC <- data.frame(Datasets,Metrics,Percentage)

# Grouped
ggplot(alignQC, aes(fill = Metrics, y = Percentage, x = Datasets)) +
  geom_bar(position = "dodge", stat = "identity", color = "black") +  # Add a black outline to the bars
  geom_text(aes(label = Percentage), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  scale_fill_manual(values = c("#dedede", "#dedede")) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(size = 12, color = "black"),  # Customize x-axis text
    axis.text.y = element_text(size = 12, color = "black"),  # Customize y-axis text
    axis.title.x = element_text(size = 14, color = "black"),  # Customize x-axis title
    axis.title.y = element_text(size = 14, color = "black")   # Customize y-axis title
  ) +
  labs(
    x = "X-Axis Label",  # Customize the x-axis label
    y = "Y-Axis Label")+  # Customize the y-axis label) +
  coord_flip()


##### Alpha Plot #####
Metrics <- rep(c("Alignment \n    Rate" , "Percentage \n  Whitelist" ), )
Percentage <- c(63.4, 95.54)

alignalphaQC <- data.frame(Metrics,Percentage)


ggplot(alignalphaQC, aes(y = Percentage, x = Metrics)) +
  geom_bar(position = "dodge", stat = "identity", fill = "#d3d3d3", color = "white") +  # Add black outline and change fill color
  geom_text(aes(label = Percentage), 
            position = position_dodge(width = 0.9), 
            hjust = -0.5,  # Adjust the horizontal position to the right of the bars
            vjust = 0.5,   # Center the labels vertically within the bars
            size = 4,  # Adjust the font size for the labels
            color = "black") +  # Adjust the font color for the labels
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove y-axis major grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.text.x = element_text(size = 12, color = "black"),  # Customize x-axis text
    axis.text.y = element_text(size = 12, color = "black", hjust = 0),  # Customize y-axis text and adjust alignment
    axis.title.x = element_text(size = 14, color = "black"),  # Customize x-axis title
    axis.title.y = element_text(size = 14, color = "black")   # Customize y-axis title
  ) +
  labs(
    x = "",  # Customize the x-axis label
    y = "") +
  coord_flip()

##### Lambda Plot #####


##### Untreated Plot ####

