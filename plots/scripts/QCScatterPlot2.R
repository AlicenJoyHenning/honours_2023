# QC Scatter Plot 


treatments <- c(rep("alpha", 4), rep("lambda", 4), rep("untreated", 4))
QCMetric <- rep(c( "Percent MT", "HIGH nF","Accepted"), 3)
cells <- c(0, 330, 1072, 4767, 0, 476, 947, 4511,0, 399, 1430, 4656)
grouped <- data.frame(treatments, QCMetric, cells)
labels <- c("0%", "5%", "18%", "77%","0%", "8%", "16%", "76%","0%", "6%", "22%", "72%")
colours <- c("LOW nF" = "white", "Percent Mt"="lightgrey", "HIGH nF" ="grey", "Accepted" = "#6ab5ba")


horizontal <-   
  ggplot(
    grouped,
    aes(fill=QCMetric, y=treatments, x=cells)) + 
  geom_bar(color = NA, position="stack", stat="identity") + 
  geom_text(inherit.aes = TRUE, aes(label= labels), vjust=0) +
  #geom_label(label = labels, aes(fill = QCMetric), colour = "black", fontface = "bold",  label.padding = unit(0, "mm")) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = colours)+ 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12, face = "bold"),   # Adjust x-axis text size and style
    axis.text.y = element_text(size = 12),                 # Adjust y-axis text size
    axis.title = element_text(size = 14, face = "bold"),   # Adjust axis title size and style
    legend.title = element_text(size = 14, face = "bold"),                        # Remove legend title
    legend.text = element_text(size = 12))                 # Adjust legend text size

