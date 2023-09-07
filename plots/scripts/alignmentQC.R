# Alignment Quality 

library(ggplot2)
install.packages("extrafont")
library(extrafont)
font_import(pattern = "EBGaramond")

Datasets <- c(rep("alpha" , 2), rep("lambda", 2), rep("untreated", 2)) 
Metrics <- rep(c("Alignment Rate" , "Percentage Whitelist" ),3)
Percentage <- c(56.6, 95.6, 55.3, 95.7, 54.1, 96.1)

alignQC <- data.frame(Datasets,Metrics,Percentage)

# Grouped
ggplot(alignQC, aes(fill=Metrics, y=Percentage, x=Datasets)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label = Percentage),  position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("#dedede", "#6ab5ba")) +
 # NO LITERATURE TO SUPPORT geom_hline(yintercept =50, linetype = "dotted", color = "black") + 
 # geom_hline(yintercept = 98, linetype = "dotted", color = "darkblue") +
  coord_flip()


ggplot(data, aes(fill = Metrics, y = Percentage, x = Datasets)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("#dedede", "#6ab5ba")) +
  coord_flip() +
  
  # Add labels inside the bars at the top
  geom_text(aes(label = Percentage, y = Percentage + 5),  # Adjust the 'y' position
            position = position_dodge(width = 0.9), 
            size = 3.5, hjust = 0.5) +  # Center the labels horizontally
  
  labs(y = "Percentage") +
  theme(text = element_text(family = "EBGaramond"))


