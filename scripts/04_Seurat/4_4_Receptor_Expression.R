# Receptor expression across cell types 
setwd("..")
getwd()
# [1] Load integrated data set :
TreatmentAnnotated <- read_rds("alice/honours/results/FinalIndex/TreatmentAnnotated.rds")
TreatmentAnnotated@meta.data$cell_type <- TreatmentAnnotated@active.ident
Idents(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated) <- "RNA"

TreatmentAnnotated 

# [2] View receptor expression across cell types : 

levels(TreatmentAnnotated)
TreatmentAnnotatedPlot <- RenameIdents(TreatmentAnnotated,
                                   'monocytes' = 'mono',
                                   'naive CD4 T' = 'nCD4 T',
                                   'neutrophils'= 'neu',
                                   'T helper' = 'Th',
                                   'naive CD8 T'= 'nCD8 T',
                                   'cytotoxic T'='cCD8 T')

TreatmentAnnotated <- subset(TreatmentAnnotated, seurat_clusters != 13)

# Modify the order of CellTypes as a factor: (prevents alphabetically losing NB information)
TreatmentAnnotatedPlot$cell_type <- factor(TreatmentAnnotatedPlot$cell_type, levels = c(
  "monos", "neu", "DCs","mDCs","pDCs","platelets",
  "Th","nCD4 T","nCD8 T","cCD8 T","NKT","Tregs","Tcm","NK",
  "B")) 

levels(TreatmentAnnotatedPlot) <- #c(0,4,2,14,8,17,12,3,1,5,7,9,10,15,16,11,6,13)
  c("mono", "neu", "DCs","mDCs","pDCs","platelets",
   "Th","nCD4 T","nCD8 T","cCD8 T","NKT","Tregs","Tcm","NK",
   "B", "unknown")

AreceptorDP <- DotPlot(TreatmentAnnotated,
                      features = c("IFNAR1"), # ,  "IFNLR1"
                      cols = c("#a9aaa9", "#a9aaa9")) +
  labs(x = "", y = "", fill = "Gene expression (%)") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = c(4.65, 5.35, 14.65, 15.35), color = "black", linetype = "dashed", size = 0.25)# +

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

ReceptorTreatment <- subset(ReceptorTreatment, seurat_clusters != 13)

receptorDP <- DotPlot(TreatmentAnnotatedPlot,
                      features = c("IFNAR1","IFNAR2", "IFNLR1", "IFNLR2"),
                      cols = c("black", "black")) + ##a9aaa9
  labs(x = "", y = "") +
  theme(
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom") +
  geom_hline(yintercept = c(4.75, 5.25, 14.75, 15.25), color = "black", linetype = "dashed", size = 0.25) +
  coord_flip()

print(receptorDP)
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
  