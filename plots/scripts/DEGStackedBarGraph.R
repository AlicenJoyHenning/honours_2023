# DEG quantification stacked bar graph 

##### Stacked Bar graph attempt #####
# Y values : 
DEGs <- c(1360, 67, 
          1413, 286, 
          1507, 175, 
          762, 28, 
          951, 37, 
          1465, 602, 
          961, 69, 
          2160, 923, 
          1287, 310)
# x groups : 
CellTypes <- c(rep("cMono", 2), 
               rep("intMono", 2),
               rep("neu", 2),
               rep("CD4h", 2),
               rep("CD4n", 2),
               rep("Tregs", 2),
               rep("CD8", 2),
               rep("NK", 2),
               rep("B", 2))
# stacks : 
stimulation <- rep(c("alpha", "lambda"), 9)
# other : 
grouped <- data.frame(stimulation, CellTypes, DEGs)

colours <- c("alpha" = "#c35cad", "lambda" = "#e1aad4",
             "alpha" = "#6ab5ba", "lambda" ="#b0d8da",
             "alpha" = "#2e8f95","lambda" ="#86dbe0",
             "alpha" = "#8caf2e","lambda" ="#bfda71",
             "alpha" = "#ee5e17","lambda" ="#f49c6e",
             "alpha" = "#d72554","lambda" ="#eb8da7",
             "alpha" =  "#900c3e","lambda" ="#ac4a70",
             "alpha" =  "#265221","lambda" ="#61b061",
             "alpha" = "#00a68e","lambda" ="#71cdc0"
)

# colours <- c(
#   "cMono" = c(alpha = "#c35cad", lambda = "#e1aad4"),
#   "intMono" = c(alpha = "#6ab5ba", lambda = "#b0d8da"),
#   "neu" = c(alpha = "#2e8f95", lambda = "#86dbe0"),
#   "CD4h" = c(alpha = "#8caf2e", lambda = "#bfda71"),
#   "CD4n" = c(alpha = "#ee5e17", lambda = "#f49c6e"),
#   "Tregs" = c(alpha = "#d72554", lambda = "#eb8da7"),
#   "CD8" = c(alpha = "#900c3e", lambda = "#ac4a70"),
#   "NK" = c(alpha = "#265221", lambda = "#61b061"),
#   "B" = c(alpha = "#00a68e", lambda = "#71cdc0")
# )

# Plot : 
horizontal <- ggplot(
  grouped,
  aes(fill = stimulation, y = DEGs, x = CellTypes)) +
  geom_bar(color = "black", position = "stack", stat = "identity") +
  geom_label(label = DEGs, aes(fill = stimulation), colour = "black", fontface = "bold", size = 5.5) +
  theme_minimal() +
  scale_fill_manual(values = colours) +
  labs(title = "", x = "Immune Cell Types", y = "DEGs") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  )


##### Giving up and making it a normal bar graph #####
# Y values : 
DEGs <- c(1360, 67, 
          1413, 286, 
          1507, 175, 
          762, 28, 
          951, 37, 
          1465, 602, 
          961, 69, 
          2160, 923, 
          1287, 310)
# x groups : 
CellTypes <- c("A_cMono","L_cMono", 
               "A_intMono","L_intMono",
               "A_neu", "L_neu",
               "A_CD4h","L_CD4h",
               "A_CD4n", "L_CD4n",
               "A_Tregs", "L_Tregs",
               "A_CD8", "L_CD8",
               "A_NK", "L_NK",
               "A_B","L_B")

colours <- c("#c35cad", "#e1aad4",
             "#6ab5ba", "#b0d8da",
             "#2e8f95","#86dbe0",
             "#8caf2e","#bfda71",
             "#ee5e17","#f49c6e",
             "#d72554","#eb8da7",
             "#900c3e","#ac4a70",
             "#265221","#61b061",
             "#00a68e","#71cdc0")

# Create a data frame
grouped <- data.frame(CellTypes, DEGs)


# Plot
horizontal <- ggplot(
  grouped,
  aes(x = CellTypes, y = DEGs, fill = CellTypes)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colours) +
  theme_minimal() +
  labs(title = "", x = "Immune Cell Types", y = "DEGs") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none"  # Hide the legend
  )

horizontal



