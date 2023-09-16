# DEG quantification stacked bar graph 

##### Stacked Bar graph attempt #####
# Y values : 
DEGs <- c(442, 17, -726, -7,
          316, 6, -199, -3, 
          418, 19, -305, -7, 
          438, 26, -357, -7,
          126, 6, -12, -2,
          472, 44, -382, -7,
          194, 2, -76, -1,
          337, 79, -179, -3
)
# x groups : 
CellTypes <- c(rep("Mono", 4), 
               rep("Neu", 4),
               rep("CD4h", 4),
               rep("CD4n", 4),
               rep("Tregs", 4),
               rep("CD8", 4),
               rep("NK", 4),
               rep("B", 4))
# stacks : 
stimulation <- rep(c("alpha upreg", "lambda upreg", "alpha downreg", "lambda downreg"), 8)
# other : 
grouped <- data.frame(stimulation, CellTypes, DEGs)

colours <- c("#b0d8da", "#57a8eb","#6ab5ba","#187bcd")
colours <- c("#d3d3d3","#d3d3d3", "#6ab5ba","#6ab5ba")

# Plot : 
horizontal <- ggplot(
  grouped,
  aes(fill = stimulation, y = DEGs, x = CellTypes)) +
  geom_bar(color = "white", position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = colours) +
  labs(title = "", x = "Immune Cell Types", y = "DEGs") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

vertical <- ggplot(
  grouped,
  aes(fill = stimulation, y = CellTypes, x = DEGs)) +
  geom_bar(color = "white", position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = colours) +
  labs(title = "", x = "DEGs", y = "Immune Cell Types") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


##### just DEs, not up or down #####

# Y values : 
DEGs <- c(1168, 24,
          515, 9,
          723, 26, 
          795, 33, 
          138, 8, 
          854, 51, 
          270, 3, 
          516, 82
)
# x groups : 
CellTypes <- c(rep("Mono", 2), 
               rep("Neu", 2),
               rep("CD4h", 2),
               rep("CD4n", 2),
               rep("Tregs", 2),
               rep("CD8", 2),
               rep("NK", 2),
               rep("B", 2))
# stacks : 
Treatment <- rep(c("alpha", "lambda"), 8)
# other : 
grouped <- data.frame(Treatment, CellTypes, DEGs)

colours <- c("lightgrey", "#6ab5ba")

# Plot : 
horizontal <- ggplot(
  grouped,
  aes(fill = Treatment, y =DEGs, x = CellTypes)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = colours) +
  labs(title = "", x = "Number of DEGs", y = "Immune Cell Types") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


##### DGEA #####

# L6: B cells lambda specific 
# Y values : 
process <- L6_lambda$Description
# x groups : 
Count <- L6_lambda$Count

pvalue <- L6_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
B_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "B", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot")



# L2: CD4 naive lambda specific 
# Y values : 
process <- L2_lambda$Description
# x groups : 
Count <- L2_lambda$Count

pvalue <- L2_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
L2_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "CD4 naive", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot")


# L1: CD4 helper lambda specific 
# Y values : 
process <- L1_lambda$Description
# x groups : 
Count <- L1_lambda$Count

pvalue <- L1_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
L1_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "CD4 helper", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot")

# M2: neutrophils lambda specific 

subset <- c("GO:0006775", "GO:0071354","GO:0033280","GO:0071295","GO:0071305","GO:0042359","GO:0009111","GO:0042363","GO:1990573","GO:0070741")

M2_lambda <- M2_lambda[M2_lambda$ID %in% subset,]

# Y values : 
process <- M2_lambda$Description
# x groups : 
Count <- M2_lambda$Count

pvalue <- M2_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
M2_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "Neutrophils", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot")

# L4: CD8 lambda specific
L4_lambda$Description[5] <- "DNA damage response signalling"

# Y values : 
process <- L4_lambda$Description
# x groups : 
Count <- L4_lambda$Count

pvalue <- L4_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
L4_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack", width = 1) + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "CD8", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot") 


# L3  Tregs lambda specific 

subset <- c("GO:0002181", "GO:0042273", "GO:0030643", "GO:0042363", "GO:0009111","GO:0044320", "GO:0042359","GO:0044321","GO:0071305","GO:0006706","GO:0033280", "GO:0030279", "GO:0071354", "GO:0006775", "GO:0031670","GO:0070741", "GO:0055062", "GO:0071107", "GO:0044320")

L3_lambda <- L3_lambda[L3_lambda$ID %in% subset,]

# Y values : 
process <- L3_lambda$Description
# x groups : 
Count <- L3_lambda$Count

pvalue <- L3_lambda$p.adjust

grouped <- data.frame(process, Count, pvalue)

# Plot
L3_bar <- ggplot(
  grouped,
  aes(fill = pvalue, y = process, x = Count)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + 
  theme_minimal() +
  scale_fill_gradient(low = "#d3d3d3", high = "#6ab5ba") +  # Use scale_fill_gradient
  labs(title = "Tregulatory", x = "Count", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    title = element_text(size = 20, face = "bold", colour = "black", vjust = 0.5),
    axis.title = element_text(size = 18, face = "bold", colour = "black", margin = margin(t = 10)),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title.position = "plot")

