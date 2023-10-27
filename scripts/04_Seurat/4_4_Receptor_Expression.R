# Receptor expression across cell types 
library(readr)
library(Seurat)
library(ggplot2)

BiocManager::install("dynverse")


setwd("..")
getwd()
# [1] Load integrated data set :
treatment <- readRDS("honours/results/FinalIndex/adjtreatment.rds")
Idents(treatment)
DefaultAssay(treatment) <- "RNA"
plot <- RenameIdents(treatment,
                     '0' = 'monocytes',
                                   '1' = 'naive CD4 T',
                                   '2' = 'neutrophils',
                                   '3' = 'T helper',
                                   '4' = 'monocytes',
                                   '5' = 'naive CD8 T',
                                   '6' = 'B',
                                   '7' = 'cytotoxic T',
                                   '8' = 'mDCs',
                                   '9' = 'NKT',
                                   '10' = 'Tregs',
                                   '11' = 'NK',
                                   '12' = 'platelets',
                                   '13' = 'unknown',
                                   '14' = 'DCs', 
                                   '15' = 'Tcm',
                                   '16' = 'Tcm',
                                   '17' = 'pDCs')
plot@meta.data$celltype <- plot@active.ident
plot <- subset(plot, celltype != "platelets" & celltype != "unknown")


# [2] View receptor expression across cell types : 

# ABBREVIATED 
# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
# TreatmentAnnotated$celltype <- factor(TreatmentAnnotated$celltype, levels = c(
#   "monos", "neu", "DCs","mDCs","pDCs",
#   "Th","nCD4 T","nCD8 T","cCD8 T","NKT","Tregs","Tcm","NK",
#   "B")) 
# 
# levels(TreatmentAnnotated) <- #c(0,4,2,14,8,17,12,3,1,5,7,9,10,15,16,11,6,13)
#   c("mono", "neu", "DCs","mDCs","pDCs",
#    "Th","nCD4 T","nCD8 T","cCD8 T","NKT","Tregs","Tcm","NK",
#    "B")

# FULL NAMES
# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
plot$celltype <- factor(plot$celltype, levels = c(
  "monocytes", "neutrophils", "DCs", "mDCs","pDCs", # "platelets",
  "T helper","naive CD4 T","naive CD8 T","cytotoxic T","NKT","Tregs","Tcm","NK",
  "B"))
levels(plot) <- #c(0,4,2,14,8,17,12,3,1,5,7,9,10,15,16,11,6,13)
  c( "monocytes", "neutrophils", "DCs", "mDCs","pDCs",# "platelets",
     "T helper","naive CD4 T","naive CD8 T","cytotoxic T","NKT","Tregs","Tcm","NK",
"B")

# PLOT
DotPlot(plot,
        features = c("IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1", "IL10RB"), # ,  
        cols = c("darkgrey", "darkgrey")) +
  labs(x = "", y = "", fill = "Gene expression (%)") +
  theme(
    axis.text.x = element_text( hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = c(4.65, 5.35), color = "black", linetype = "dashed", size = 0.25)# +

LreceptorDP <- DotPlot(TreatmentAnnotated,
                       features = c("IFNLR1"), # ,  "IFNLR1"
                       cols = c("#6ab5ab", "#6ab5ab")) +
  labs(x = "", y = "", fill = "Gene Expression (%)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    legend.title = element_text(size = 14, colour= "black", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  geom_hline(yintercept = c(4.65, 5.35, 14.65, 15.35), color = "black", linetype = "dashed", size = 0.25)# +

AreceptorDP /   LreceptorDP                        
#blend to view alpha and lambda receptors!
?FeaturePlot()

# [3] Combining T cell populations 

ReceptorTreatment <- RenameIdents(TreatmentAnnotated, 
                                  'naive CD4 T' = "T cells",
                                  'T helper' = "T cells",
                                  'naive CD8 T' = "T cells",
                                  'cytotoxic T' = "T cells", 
                                  'NKT' = "T cells",
                                  'Tregs' = "T cells",
                                  'NK' = "T cells",
                                  'Tcm' = "T cells",
                                  'Tcm'= "T cells",
                                  "B" = 'B cells')

ReceptorTreatment <- subset(ReceptorTreatment, celltype != 13)

Breceptors <- DotPlot(TreatmentAnnotatedPlot,
                      features = c("CD79A", "CD79B"), # "IFNAR1",IFNAR2", "IFNLR1", "IL10RB"
                      cols = c("white", "grey")) + ##a9aaa9
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom") +
  geom_hline(yintercept = c(14.75, 15.25), color = "black", linetype = "dashed", size = 0.25) +
  coord_flip() # 4.75, 5.25, 

Breceptors <- DotPlot(Bcells,
                      features = c("ISGF3"),
                      split.by = "treatment",
                      cols = c("grey", "grey", "black")) + ##a9aaa9
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom") +
  geom_hline(yintercept = c(4.75, 5.25, 14.75, 15.25), color = "black", linetype = "dashed", size = 0.25)+
coord_flip() +
  scale_y_discrete(labels = c("u", "λ", "α"))

IFNreceptors / Breceptors

# Plot : 

receptorexpression <- ggplot(
  grouped,
  aes(fill = rev(Treatment), y = CellTypes, x = Goterms)) +
  geom_bar(color = NA, position = "stack", stat = "identity") +
  theme_minimal() +
  scale_fill_manual(values = rev(colours)) +
  labs(title = "", x = "Number of cells", y = "", fill = "IFN treatment") +
  theme(
    axis.text.x = element_text(size = 16, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black", hjust = 0.5),  # Center-align y-axis text
    axis.title.x =  element_text(size = 16, face = "bold", colour = "black", margin = margin(t = 10)),
    axis.title.y =  element_blank(),legend.title = element_text(size = 14, face = "bold", colour= "black", margin = margin(t = -50)),
    legend.text = element_text(size = 14, colour= "black"),
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

citation(package = "Seurat")
