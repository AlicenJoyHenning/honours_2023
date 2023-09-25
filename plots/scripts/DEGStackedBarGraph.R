# DEG quantification stacked bar graph 

library(ggplot2)

##### Stacked Bar graph DEGs up & down #####
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
  #labs(title = "", x = "DEGs", y = "Immune Cell Types") +
  labs(title = "", x = "DEGs", y = "enrichGO terms") +
  theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


##### Stacked Bar graph DEGs #####

# Y values : 
DEGs <- c(1028, 10, 1, # monocytes
          699, 4, 6, # neutrophils
          # 1209, 13, # myeloid 
          3, 0, 0,   # DCs
          80, 0, 0, # myeloid dericed dendritic cells
          0, 0, 0, # pDCs

          # 297, 51, # dendritic cells 

          12, 0, 0, # platelets 

          807, 12, 11, # T cells
          714, 13, 13, # CD4h
          678, 5, 6, # naive cd8
          310, 5, 3, # NKT
          186, 1, 1, # cyto CD8
          148, 5, 1, # T regs 
          62, 0, 0, # CD4
          150, 1, 1, # NK
          # 935, 33, # overall T 

          350, 2, 34 # B
)

# x groups : 
CellTypes <- c(rep("mono", 3), 
               rep("neu", 3),
               # rep("myeloid", 2),
               rep("DCs", 3),
               rep("mDCs", 3),
               rep("pDCs", 3),

               rep("platelets", 3),

               rep("T", 3),
               rep("CD4h", 3),
               rep("nCD8", 3),
               rep("NKT", 3),
               rep("cCD8", 3),
               rep("Tregs", 3),
               rep("CD4", 3),
               rep("NK", 3),
               # rep("all_T", 2),

               rep("B", 3))
# stacks : 
Treatment <- rep(c("alpha","common","lambda"), 15)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, DEGs)

# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "mono", "neu", "DCs","mDCs","pDCs","platelets","T","CD4h","nCD8","NKT","cCD8","Tregs","CD4","NK","B"))
grouped$Treatment <- factor(grouped$Treatment, levels = c("alpha", "common", "lambda"))

colours <- c("lightgrey","#6ab5ba","darkgrey")


# Plot : 
verticalDEGs <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = -DEGs)) +  # Reverse x-axis and fill order
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +  # Reverse fill colors
  labs(title = "", x = "Number of DEGs", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank(),
    #legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  guides(fill = "none")  # Remove the legend

##### Stacked Bar graph for GO terms #####

# Y values : 
Goterms <- c(766, 0, 7, # monocytes
          598, 2, 67, # neutrophils
          # 1209, 13, # myeloid 
          53, 0, 0,   # DCs
          125, 0, 0, # myeloid dericed dendritic cells
          0, 0, 0, # pDCs
          
          # 297, 51, # dendritic cells 
          
          11, 0, 0, # platelets 
          
          493, 1, 2, # T cells
          406, 3, 6, # CD4h
          416, 4, 0, # naive cd8
          282, 1, 0, # NKT
          190, 1, 0, # cyto CD8
          203, 0, 13, # T regs 
          129, 0, 0, # CD4
          207, 0, 1, # NK
          # 935, 33, # overall T 
          
          339, 70, 10 # B
)

# x groups : 
CellTypes <- c(rep("mono", 3), 
               rep("neu", 3),
               # rep("myeloid", 2),
               rep("DCs", 3),
               rep("mDCs", 3),
               rep("pDCs", 3),
               
               rep("platelets", 3),
               
               rep("T", 3),
               rep("CD4h", 3),
               rep("nCD8", 3),
               rep("NKT", 3),
               rep("cCD8", 3),
               rep("Tregs", 3),
               rep("CD4", 3),
               rep("NK", 3),
               # rep("all_T", 2),
               
               rep("B", 3))
# stacks : 
Treatment <- rep(c("alpha", "common", "lambda"), 15)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, Goterms)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, -Goterms)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "mono", "neu", "DCs","mDCs","pDCs","platelets","T","CD4h","nCD8","NKT","cCD8","Tregs","CD4","NK","B"))


colours <- c("lightgrey", "darkgrey", "#6ab5ba")


# Plot : 
verticalGO <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = Goterms)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +
  # number of enriched GO terms 
  labs(title = "", x = "Number of enriched GO terms", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

##### Stacked bar graph together #####
verticalDEGs | verticalGO 

##### UP & DOWN regulated myeloid #####
ALPHA
# Y values : 
alphaDEGs <- c(414, 615, # monocytes
             317, 386, # neutrophils
             2, 1,   # DCs
             59, 21, # myeloid dericed dendritic cells
             0, 0) # pDCs)

# x groups : 
CellTypes <- c(rep("mono", 2), 
               rep("neu", 2),
               rep("DCs", 2),
               rep("mDCs", 2),
               rep("pDCs", 2))

# stacks : 
Treatment <- rep(c("up regulated", "down regulated"), 5)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, alphaDEGs)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, alphaDEGs)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "mono", "neu", "DCs","mDCs","pDCs"))


colours <- c("lightgrey", "#6ab5ba")


# Plot : 
myeloidalpha <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = -alphaDEGs)) +  # Reverse x-axis and fill order
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +  # Reverse fill colors
  labs(title = "", x = "Number of DEGs (α)", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank(),
    #legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  guides(fill = "none")

# LAMBDA 
# Y values : 
lambdaDEGs <- c(10, 1, # monocytes
               9, 1, # neutrophils
               0, 0,   # DCs
               0, 0, # myeloid dericed dendritic cells
               0, 0) # pDCs)

# x groups : 
CellTypes <- c(rep("mono", 2), 
               rep("neu", 2),
               rep("DCs", 2),
               rep("mDCs", 2),
               rep("pDCs", 2))

# stacks : 
Treatment <- rep(c("up regulated", "down regulated"), 5)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, lambdaDEGs)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, -lambdaDEGs)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "mono", "neu", "DCs","mDCs","pDCs"))


colours <- c("lightgrey","#6ab5ba")


# Plot : 
myeloidlambda <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = lambdaDEGs)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +
  # number of enriched GO terms 
  labs(title = "", x = "Number of DEGs (λ)", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


myeloidalpha | myeloidlambda

##### UP & DOWN regulated lymphoid #####
ALPHA
# Y values : 
alphaDEGs <- c(504, 315, # T cells
          422, 305, # CD4h
          433, 305, # naive cd8
          204, 111, # NKT
          136, 51, # cyto CD8
          144, 9, # T regs 
          58, 4, # CD4
          131, 20, # NK
          
          567, 368 ) # B

# x groups : 
CellTypes <- c(rep("T", 2), 
               rep("CD4h", 2),
               rep("nCD8", 2),
               rep("NKT", 2),
               rep("cCD8", 2),
               rep("Tregs", 2),
               rep("CD4", 2), 
               rep("NK", 2),
               rep("B", 2))
# stacks : 
Treatment <- rep(c("up regulated", "down regulated"), 9)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, alphaDEGs)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, alphaDEGs)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "T","CD4h","nCD8","NKT", "cCD8","Tregs","CD4","NK","B"))


colours <- c("lightgrey", "#6ab5ba")


# Plot : 
lymphoidalpha <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = -alphaDEGs)) +  # Reverse x-axis and fill order
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +  # Reverse fill colors
  labs(title = "", x = "Number of DEGs (α)", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank(),
    #legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  guides(fill = "none")

# LAMBDA 
# Y values : 
lambdaDEGs <- c(22, 1, # T cells
               22, 4, # CD4h
               10, 1, # naive cd8
               6, 2, # NKT
               2, 0, # cyto CD8
               4, 2, # T regs 
               0, 0, # CD4
               1, 0, # NK
               
               30, 3 ) # B

# x groups : 
CellTypes <- c(rep("T", 2), 
               rep("CD4h", 2),
               rep("nCD8", 2),
               rep("NKT", 2),
               rep("cCD8", 2),
               rep("Tregs", 2),
               rep("CD4", 2), 
               rep("NK", 2),
               rep("B", 2))
# stacks : 
Treatment <- rep(c("up regulated", "down regulated"), 9)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, lambdaDEGs)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, lambdaDEGs)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "T","CD4h","nCD8","NKT", "cCD8","Tregs","CD4","NK","B"))

colours <- c("lightgrey","#6ab5ba")


# Plot : 
lymphoilambda <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = lambdaDEGs)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +
  # number of enriched GO terms 
  labs(title = "", x = "Number of DEGs (λ)", y = "") +
  theme(
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_text(size = 16, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 16, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )


lymphoidalpha | lymphoilambda

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

