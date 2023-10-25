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
          62, 0, # CD4
          150, 2, # NK
          # 935, 33, # overall T 
          
          350, 36 # B
)

# x groups : 
CellTypes <- c(rep("mono", 2), 
               rep("neu", 2),
               # rep("myeloid", 2),
               rep("DCs", 2),
               rep("mDCs", 2),
               rep("pDCs", 2),
               
               # rep("platelets", 3),
               
               rep("Th", 2),
               rep("nCD4+ T", 2),
               rep("nCD8+ T", 2),
               rep("cyto T", 2),
               rep("NKT", 2),
               rep("Treg", 2),
               rep("Tcm", 2),
               rep("NK", 2),
               # rep("all_T", 2),
               
               rep("B", 2))
# stacks : 
Treatment <- rep(c("alpha","lambda"), 14)

# create data frame :  
grouped <- data.frame(Treatment, CellTypes, DEGs)

# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
grouped$CellTypes <- factor(grouped$CellTypes, levels = c(
  "mono", "neu", "DCs","mDCs","pDCs",
  "Th","nCD4+ T","nCD8+ T","cyto T","NKT","Treg","Tcm","",
  "B"))
grouped$Treatment <- factor(grouped$Treatment, levels = c("alpha", "common", "lambda"))

colours <- c("#6ab5ba","#efefef")


verticalDEGs <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = DEGs)) +
  geom_bar(color = "black", position = "stack", stat = "identity", size = 0.5) +  # Set size for the outline
  theme_minimal() +
  scale_fill_manual(values = colours) +
  # number of enriched GO terms 
  labs(title = "", x = "Number of DEGs", y = "", fill = "IFN treatment") +
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
