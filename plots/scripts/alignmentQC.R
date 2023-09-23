# Alignment Quality 

##### Loading dependencies ####
library(ggplot2)
install.packages("extrafont")
library(extrafont)
font_import(pattern = "EBGaramond")

# BLUE "#a0d9e9"

##### Merged Plot (all datasets) ####
Datasets <- c(rep("α" , 2), rep("λ", 2), rep("u", 2)) 
Metrics <- rep(c("Percentage Psuedoaligned" , "Percentage Whitelist" ),3)
Percentage <- c(63.4, 95.5, 61.4, 95.6, 60.5, 96.0)

alignQC <- data.frame(Datasets,Metrics,Percentage)

# Grouped
ggplot(alignQC, aes(fill = Metrics, y = Percentage, x = Datasets)) +
  geom_bar(position = "dodge", stat = "identity", color = "white") +  # Add a black outline to the bars
  geom_text(aes(label = Percentage), position = position_dodge(width = 0.9), vjust = -0.5, hjust = 0.5, size = 5.5) +
  theme_minimal() +
  scale_fill_manual(values = c("#dedede", "grey")) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(size = 18,  color = "black"),  # Adjust x-axis text angle and hjust ),  # Customize x-axis text
    axis.text.y = element_text(size = 18, color = "black"),  # Customize y-axis text
    axis.title.x = element_text(size = 20, color = "black", vjust = -0.2),  # Customize x-axis title
    axis.title.y = element_text(size = 20, color = "black",  face = "bold"),  # Customize y-axis title
    legend.text = element_text(size = 18),  # Adjust legend text size
    legend.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "",  # Customize the x-axis label
    y = "Percentage",  # Customize the y-axis label
    fill = "Barcode Metrics") + # customize legend 
  scale_x_discrete(expand = c(0.5, 0))   # Adjust x-axis limits to expand the range
  #coord_flip() 


##### Alpha Plot #####
Metrics <- rep(c("Alignment \n    Rate" , "Percentage \n  Whitelist" ), )
Percentage <- c(63.4, 95.54)

alignalphaQC <- data.frame(Metrics,Percentage)


aligna <- 
  ggplot(alignalphaQC, aes(y = Percentage, x = Metrics)) +
  geom_bar(position = "dodge", stat = "identity", fill = "#d3d3d3", color = "white") +  # Add black outline and change fill color
  geom_text(aes(label = Percentage), 
            #position = position_dodge(width = 0.9), 
            hjust = -0.5,  # Adjust the horizontal position to the right of the bars
            #vjust = -.5,   # Center the labels vertically within the bars
            size = 8,  # Adjust the font size for the labels
            color = "black") +  # Adjust the font color for the labels
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove y-axis major grid lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.text.x = element_text(size = 20, color = "black"),  # Customize x-axis text
    axis.text.y = element_text(size = 20, color = "black", hjust = 0),  # Customize y-axis text and adjust alignment
    axis.title.x = element_text(size = 20, color = "black"),  # Customize x-axis title
    axis.title.y = element_text(size = 20, color = "black")   # Customize y-axis title
  ) +
  labs(
    x = "",  # Customize the x-axis label
    y = "",
    legend =) +
  coord_flip()

##### Lambda Plot #####
Metrics <- rep(c("Alignment \n    Rate" , "Percentage \n  Whitelist" ), )
Percentage <- c(61.4, 95.61)

lambdaQC <- data.frame(Metrics,Percentage)


alignl <- 
  ggplot(lambdaQC, aes(y = Percentage, x = Metrics)) +
  geom_bar(position = "dodge", stat = "identity", fill = "#d3d3d3", color = "white") +  # Add black outline and change fill color
  geom_text(aes(label = Percentage), 
            #position = position_dodge(width = 0.9), 
            hjust = -0.5,  # Adjust the horizontal position to the right of the bars
            #vjust = -.5,   # Center the labels vertically within the bars
            size = 8,  # Adjust the font size for the labels
            color = "black") +  # Adjust the font color for the labels
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove y-axis major grid lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.text.x = element_text(size = 20, color = "black"),  # Customize x-axis text
    axis.text.y = element_text(size = 20, color = "black", hjust = 0),  # Customize y-axis text and adjust alignment
    axis.title.x = element_text(size = 20, color = "black"),  # Customize x-axis title
    axis.title.y = element_text(size = 20, color = "black")   # Customize y-axis title
  ) +
  labs(
    x = "",  # Customize the x-axis label
    y = "") +
  coord_flip()

##### Untreated Plot ####
Metrics <- rep(c("Alignment \n    Rate" , "Percentage \n  Whitelist" ), )
Percentage <- c(60.5, 95.98)

untreatedQC <- data.frame(Metrics,Percentage)

alignu <- 
  ggplot(untreatedQC, aes(y = Percentage, x = Metrics)) +
  geom_bar(position = "dodge", stat = "identity", fill = "#d3d3d3", color = "white") +  # Add black outline and change fill color
  geom_text(aes(label = Percentage), 
            #position = position_dodge(width = 0.9), 
            hjust = -0.5,  # Adjust the horizontal position to the right of the bars
            #vjust = -.5,   # Center the labels vertically within the bars
            size = 8,  # Adjust the font size for the labels
            color = "black") +  # Adjust the font color for the labels
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove y-axis major grid lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    axis.text.x = element_text(size = 20, color = "black"),  # Customize x-axis text
    axis.text.y = element_text(size = 20, color = "black", hjust = 0),  # Customize y-axis text and adjust alignment
    axis.title.x = element_text(size = 20, color = "black"),  # Customize x-axis title
    axis.title.y = element_text(size = 20, color = "black")   # Customize y-axis title
  ) +
  labs(
    x = "",  # Customize the x-axis label
    y = "") +
  coord_flip()

# altogether : 

aligna / alignl / alignu
