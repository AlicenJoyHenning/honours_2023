# FINAL DEG quantification stacked bar graph 

library(ggplot2)

# Y values : 
DEGs <- c(1028, 11, # monocytes
          699, 10, # neutrophils
          # 1209, 13, # myeloid 
          3, 0,   # DCs
          80, 0, # myeloid dericed dendritic cells
          0, 0, # pDCs
          
          # 297, 51, # dendritic cells 
          
          # 12, 0, 0, # platelets 
          
          807, 23, # T helper
          714, 26, # naive CD4
          678, 11, # naive cd8
          310, 8, # NKT
          186, 2, # cyto CD8
          148, 6, # T regs 
          62, 0, # Tcm
          150, 2, # NK
          # 935, 33, # overall T 
          
          350, 36 # B
)

# x groups : 
CellTypes <- c(rep("monocytes", 2), 
               rep("neutrophils", 2),
               # rep("myeloid", 2),
               rep("DCs", 2),
               rep("mDCs", 2),
               rep("pDCs", 2),
               
               # rep("platelets", 3),
               
               rep("T helper", 2),
               rep("naive CD4+ T", 2),
               rep("naive CD8+ T", 2),
               rep("cytotoxic T", 2),
               rep("natural killer T", 2),
               rep("T regulatory", 2),
               rep("T central memory", 2),
               rep("natural killer", 2),
               # rep("all_T", 2),
               
               rep("B", 2))
# stacks : 
Treatment <- rep(c("alpha","lambda"), 14)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, DEGs)

# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "monocytes", "neutrophils", "DCs","mDCs","pDCs",
  "T helper","naive CD4+ T","naive CD8+ T","cytotoxic T","natural killer T","T regulatory","T central memory","natural killer",
  "B"))
grouped$Treatment <- factor(grouped$Treatment, levels = c("alpha", "common", "lambda"))

colours <- c("#efefef","#6ab5ba")


# Plot : 
verticalDEGs <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = -DEGs)) +  # Reverse x-axis and fill order
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +  # Reverse fill colors
  labs(title = "", x = "Number of DEGs", y = "") +
  theme(
    axis.text.x = element_text(size = 16, colour = "black"),
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

verticalDEGs <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = -DEGs)) +  # Reverse x-axis and fill order
  geom_bar(color = "black", position = "stack", stat = "identity", size = 0.5) +  # Set size for the outline
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +  # Reverse fill colors
  labs(title = "", x = "Number of DEGs", y = "") +
  theme(
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.title.x =  element_text(size = 24, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_blank(),
    legend.text = element_blank(),
    #legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) + 
  guides(fill = "none")




##### Stacked Bar graph for GO terms #####

# Y values : 
Goterms <- c(766, 7, # monocytes
             598, 69, # neutrophils
             # 1209, 13, # myeloid 
             53, 0,   # DCs
             125, 0, # myeloid dericed dendritic cells
             0, 0, # pDCs
             
             # 297, 51, # dendritic cells 
             
             #11, 0, 0, # platelets 
             
             493, 3, # T cells
             406, 9, # CD4h
             416, 4, # naive cd8
             282, 1, # NKT
             190, 1, # cyto CD8
             203, 13, # T regs 
             129, 0, # CD4
             207, 1, # NK
             # 935, 33, # overall T 
             
             339, 80 # B
)

# x groups : 
# x groups : 
CellTypes <- c(rep("monocytes", 2), 
               rep("neutrophils", 2),
               # rep("myeloid", 2),
               rep("DCs", 2),
               rep("mDCs", 2),
               rep("pDCs", 2),
               
               # rep("platelets", 3),
               
               rep("T helper", 2),
               rep("naive CD4+ T", 2),
               rep("naive CD8+ T", 2),
               rep("cytotoxic T", 2),
               rep("NKT", 2),
               rep("Tregs", 2),
               rep("Tcm", 2),
               rep("NK", 2),
               # rep("all_T", 2),
               
               rep("B", 2))
# stacks : 
Treatment <- rep(c("alpha", "lambda"), 14)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, Goterms)

# Reorder the data frame based on DEGs within each CellTypes group
grouped <- grouped %>% arrange(CellTypes, Treatment, -Goterms)


# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$rownames <- factor(grouped$CellTypes, levels = c(
  "monocytes", "neutrophils", "DCs","mDCs","pDCs", #"platelets",
  "T helper","naive CD4+ T","naive CD8+ T","cytotoxic T","NKT","Tregs","Tcm","NK",
  "B"))

colours <- c("#efefef","#6ab5ba")


# Plot : 
verticalGO <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = Goterms)) +
  geom_bar(color = "black", position = "stack", stat = "identity", size = 0.5) +  # Set size for the outline
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +
  # number of enriched GO terms 
  labs(title = "", x = "Number of enriched GO terms", y = "", fill = "IFN treatment") +
  theme(
    axis.text.x = element_text(size = 24, colour = "black"),
    axis.text.y = element_text(size = 26, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 24, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),
    legend.title = element_text(size = 24, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 24, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

##### Stacked bar graph together #####
verticalDEGs | verticalGO 

