# MANUAL CLUSTER ANNOTATION 
# Seurat workflow after integration

##### [1] Load dependencies datasets #####
getwd()

library(BiocManager)
BiocManager::install('limma') # BiocManager::install('limma')
BiocManager::install("SeuratData")
BiocManager::install('multtest')
BiocManager::install('metap')
BiocManager::install('xlsx')
BiocManager::install('pheatmap')
BiocManager::install('XLConnect')
BiocManager::install('writexl')


library(dplyr)
library(readr)
library(writexl)
library(openxlsx)
library(ggplot2)
library(grid)
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(patchwork)
library(readr)
library(limma)
library(metap)
library(readxl)
library(xlsx)

# To read in the saved Seurat objects : 
treatment <- readRDS("honours/work/1109/treatment.rds")

##### [2] Identify Conserved cell type markers for cell type annotation #####

# Find markers for each cluster :
DefaultAssay(treatment) <- "RNA"

ClusterPath <- "honours/results/FinalIndex/Annotation/ClusterAnnotation.xlsx"

Cluster0Markers <- FindConservedMarkers(treatment, ident.1= 0, grouping.var = "treatment", verbose = FALSE)
Cluster0Markers <- subset(Cluster0Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster0Markers <- subset(Cluster0Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C0Rows <- rownames(Cluster0Markers)
Cluster0Markers <- cbind(RowNames = C0Rows, Cluster0Markers)
dim(Cluster0Markers) # 1992 > 1549 > 149


Cluster1Markers <- FindConservedMarkers(treatment, ident.1= 1, grouping.var = "treatment", verbose = FALSE)
Cluster1Markers <- subset(Cluster1Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster1Markers <- subset(Cluster1Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C1Rows <- rownames(Cluster1Markers)
Cluster1Markers <- cbind(RowNames = C1Rows, Cluster1Markers)
dim(Cluster1Markers) # 1645 > 1146 > 96

Cluster2Markers <- FindConservedMarkers(treatment, ident.1= 2, grouping.var = "treatment", verbose = FALSE)
Cluster2Markers <- subset(Cluster2Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster2Markers <- subset(Cluster2Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C2Rows <- rownames(Cluster2Markers)
Cluster2Markers <- cbind(RowNames = C2Rows, Cluster2Markers)
dim(Cluster2Markers) # 1885 > 1385 > 74 


Cluster3Markers <- FindConservedMarkers(treatment, ident.1= 3, grouping.var = "treatment", verbose = FALSE)
Cluster3Markers <- subset(Cluster3Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster3Markers <- subset(Cluster3Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C3Rows <- rownames(Cluster3Markers)
Cluster3Markers <- cbind(RowNames = C3Rows, Cluster3Markers)
dim(Cluster3Markers) # 1464 > 962 > 11

Cluster4Markers <- FindConservedMarkers(treatment, ident.1= 4, grouping.var = "treatment", verbose = FALSE)
Cluster4Markers <- subset(Cluster4Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster4Markers <- subset(Cluster4Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C4Rows <- rownames(Cluster4Markers)
Cluster4Markers <- cbind(RowNames = C4Rows, Cluster4Markers)
dim(Cluster4Markers) # 1501 > 700 > 29 

Cluster5Markers <- FindConservedMarkers(treatment, ident.1= 5, grouping.var = "treatment", verbose = FALSE)
Cluster5Markers <- subset(Cluster5Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster5Markers <- subset(Cluster5Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C5Rows <- rownames(Cluster5Markers)
Cluster5Markers <- cbind(RowNames = C5Rows, Cluster5Markers)
dim(Cluster5Markers) # 1360 > 675 > 11

Cluster6Markers <- FindConservedMarkers(treatment, ident.1= 6, grouping.var = "treatment", verbose = FALSE)
Cluster6Markers <- subset(Cluster6Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster6Markers <- subset(Cluster6Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C6Rows <- rownames(Cluster6Markers)
Cluster6Markers <- cbind(RowNames = C6Rows, Cluster6Markers)
dim(Cluster6Markers) # 1514 > 458 > 48 


Cluster7Markers <- FindConservedMarkers(treatment, ident.1= 7, grouping.var = "treatment", verbose = FALSE)
Cluster7Markers <- subset(Cluster7Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster7Markers <- subset(Cluster7Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C7Rows <- rownames(Cluster7Markers)
Cluster7Markers <- cbind(RowNames = C7Rows, Cluster7Markers)
dim(Cluster7Markers) # 1322 > 341 > 28 

Cluster8Markers <- FindConservedMarkers(treatment, ident.1= 8, grouping.var = "treatment", verbose = FALSE)
Cluster8Markers <- subset(Cluster8Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster8Markers <- subset(Cluster8Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C8Rows <- rownames(Cluster8Markers)
Cluster8Markers <- cbind(RowNames = C8Rows, Cluster8Markers)
dim(Cluster8Markers) # 1360 > 100 > 42 

Cluster9Markers <- FindConservedMarkers(treatment, ident.1= 9, grouping.var = "treatment", verbose = FALSE)
Cluster9Markers <- subset(Cluster9Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster9Markers <- subset(Cluster9Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C9Rows <- rownames(Cluster9Markers)
Cluster9Markers <- cbind(RowNames = C9Rows, Cluster9Markers)
dim(Cluster9Markers) # 1330 > 211 > 39 

Cluster10Markers <- FindConservedMarkers(treatment, ident.1= 10, grouping.var = "treatment", verbose = FALSE)
Cluster10Markers <- subset(Cluster10Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster10Markers <- subset(Cluster10Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C10Rows <- rownames(Cluster10Markers)
Cluster10Markers <- cbind(RowNames = C10Rows, Cluster10Markers)
dim(Cluster10Markers) # 1196 > 307 > 6 

Cluster11Markers <- FindConservedMarkers(treatment, ident.1= 11, grouping.var = "treatment", verbose = TRUE)
Cluster11Markers <- subset(Cluster11Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster11Markers <- subset(Cluster11Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C11Rows <- rownames(Cluster11Markers)
Cluster11Markers <- cbind(RowNames = C11Rows, Cluster11Markers)
dim(Cluster11Markers) # 1437 > 265 > 73 

Cluster12Markers <- FindConservedMarkers(treatment, ident.1= 12, grouping.var = "treatment", verbose = FALSE)
Cluster12Markers <- subset(Cluster12Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster12Markers <- subset(Cluster12Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C12Rows <- rownames(Cluster12Markers)
Cluster12Markers <- cbind(RowNames = C12Rows, Cluster12Markers)
dim(Cluster12Markers) # 1372 > 207 > 99

Cluster13Markers <- FindConservedMarkers(treatment, ident.1= 13, grouping.var = "treatment", verbose = FALSE)
Cluster13Markers <- subset(Cluster13Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster13Markers <- subset(Cluster13Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C13Rows <- rownames(Cluster13Markers)
Cluster13Markers <- cbind(RowNames = C13Rows, Cluster13Markers)
dim(Cluster13Markers) # 1794 > 169 > 107 

Cluster14Markers <- FindConservedMarkers(treatment, ident.1= 14, grouping.var = "treatment", verbose = FALSE)
Cluster14Markers <- subset(Cluster14Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster14Markers <- subset(Cluster14Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C14Rows <- rownames(Cluster14Markers)
Cluster14Markers <- cbind(RowNames = C14Rows, Cluster14Markers)
dim(Cluster14Markers) # 1484 > 15 > 14

Cluster15Markers <- FindConservedMarkers(treatment, ident.1= 15, grouping.var = "treatment", verbose = FALSE)
Cluster15Markers <- subset(Cluster15Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster15Markers <- subset(Cluster15Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C15Rows <- rownames(Cluster15Markers)
Cluster15Markers <- cbind(RowNames = C15Rows, Cluster15Markers)
dim(Cluster15Markers) # 1155 > 26 > 6 

Cluster16Markers <- FindConservedMarkers(treatment, ident.1= 16, grouping.var = "treatment", verbose = FALSE)
Cluster16Markers <- subset(Cluster16Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster16Markers <- subset(Cluster16Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C16Rows <- rownames(Cluster16Markers)
Cluster16Markers <- cbind(RowNames = C16Rows, Cluster16Markers)
dim(Cluster16Markers) # 881 > 0

Cluster17Markers <- FindConservedMarkers(treatment, ident.1= 17, grouping.var = "treatment", verbose = FALSE)
Cluster17Markers <- subset(Cluster17Markers, subset = alpha_p_val_adj < 0.05) 
Cluster17Markers <- subset(Cluster17Markers, subset = alpha_avg_log2FC > 1) 
C17Rows <- rownames(Cluster17Markers)
Cluster17Markers <- cbind(RowNames = C17Rows, Cluster17Markers)
dim(Cluster17Markers) # 3936 > 1384 > 271 


# Saving the data to an excel sheet to annotate genes : 
clusters <- list('Cluster0' = Cluster0Markers, 'Cluster1' = Cluster1Markers,'Cluster2' = Cluster2Markers, 'Cluster3' = Cluster3Markers,'Cluster4' = Cluster4Markers, 'Cluster5' = Cluster5Markers, 'Cluster6' = Cluster6Markers, 'Cluster7' = Cluster7Markers, 'Cluster8' = Cluster8Markers, 'Cluster9' = Cluster9Markers, 'Cluster10' = Cluster10Markers, 'Cluster11' = Cluster11Markers,'Cluster12' = Cluster12Markers,'Cluster13' = Cluster13Markers,'Cluster14' = Cluster14Markers, 'Cluster15' = Cluster15Markers, "Cluster16" = Cluster16Markers, "Cluster17" = Cluster17Markers)
openxlsx::write.xlsx(clusters, file = "honours/results/FinalIndex/Annotation/SeuratFindConservedMarkers.xlsx")
openxlsx::write.xlsx(Cluster2v3Markers, file = "honours/results/IntegratedMarkers/Cluster2v3Markers.xlsx")

##### [3] View FindConservedMarkers() for each cluster #####

# CLUSTER 0 
features <- c("CCL3", "CCL4", "CCL4L2", "CD82", "CSF3R", "CXCL8", "FPR1", "IFITM3", "ITGAX", "SRGN", "TREM1")

DotPlot(object = treatment, features = "ZNF683", cols = c("grey", "#15c284")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_hline(yintercept = c(0, 9, 14, 15, 12, 17, 18), linetype = "dotted", color = "black")

 # CLUSTER 1 
features <- c("ABLIM1", "CA6", "CAMK4", "CCR7", "DGKA", "FAM117B", "LEF1", "LEF1-AS1","LRRN3", "MAL", "MYC", "NELL2", "NOG", "OXNAD1", "PASK", "PDK1", "SERINC5", "SH3YL1", "TCF7","TMEM204","TXK")
features <- c("NKG2A")
features <- c("ITGAM", "ITGAX", "IL3RA", "FCGR3A","SIGLEC3", "HLA-DRA", "HLA-DRB1")
DotPlot(object = treatment, features = features, cols = c("grey", "#d72554")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 2 
features <- c("CD68", "IFNAR1")
DotPlot(object = treatment, features = features, cols = c("grey", "#7ac745")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 3 

DotPlot(object = treatment, features = features, cols = c("grey", "#7e549f")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 4 

DotPlot(object = treatment, features = features, cols = c("grey", "#37a777")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 5 
features <- c("CD41", "CD154", "CD61","CD62P") 
DotPlot(object = treatment, features = features, cols = c("grey", "#fb836f")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 6 = B Cells
features <- c("CD79A", "MS4A1", "BLNK", "CD79B", "CD24", "CD19")
DotPlot(object = treatment, features = features, cols = c("grey", "#a0d9e9")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 7 

DotPlot(object = treatment, features = features, cols = c("grey", "#9b78c2")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 8 : dendritic cells 
features <- c("MS4A7", "CD86","CST3", "CD14", "CD80") 
DotPlot(object = treatment, features = features, cols = c("grey", "#6ab5ba")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 9
features <- c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B")
features <- c("NKG7", "CXCR6")
DotPlot(object = treatment, features = features, cols = c("grey", "#93aff5")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 10 : Treg
features <- c("FOXP3", "RTKN2","IKZF2", "TBC1D4","STAM", "FANK1")
DotPlot(object = treatment, features = features, cols = c("grey", "#c674bc")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 11 : NK cells 
features <- c("KLRD1","GNLY","KLRB1","NKG7","KLRK1","FCGR3A","GZMB","KLRC1","KLRF1","NCAM1")
DotPlot(object = treatment, features = features, cols = c("grey", "#81cfff")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 12

DotPlot(object = treatment, features = features, cols = c("grey", "#8edecf")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 13

DotPlot(object = treatment, features = features, cols = c("grey", "#69a923")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 14

DotPlot(object = treatment, features = features, cols = c("grey", "#00a68e")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 15

DotPlot(object = treatment, features = features, cols = c("grey", "#d73f3f")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

# CLUSTER 17

DotPlot(object = treatment, features = features, cols = c("grey", "#a0d9e9")) 
+  coord_flip() 
+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust x-axis labels angle
+  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")




##### [4] Tools for Manual Annotation checking literature markers  #####

#Cluster #
DefaultAssay(treatment) <- "RNA"

features = c("SRY", "DDX3Y", "ZFY", "TSPY", "PRY", "UTY", "XIST")
features = c("CD8A", "ZNF683", "CD8B", "GZMB", "GZMA", "CD3G", "CCR7","CX3CR1", "DSTN")
features = c("CD3G", "IL17R", "CCR4", "CX3CR1", "CD161", "CCR6", "Leu8")
FeaturePlot(object = treatment, 
            features = features,
            cols = c("lightgrey", "black"),
            label = TRUE,
            pt.size = 1.5, 
            blend = FALSE)  
  # theme(panel.background = element_rect(fill = "darkgrey")) 

cluster_boundaries <- c(5, 7)

DotPlot(object = treatment, 
        features = features)
        cols = c("grey", "#265221")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Adjust x-axis labels angle
  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

##### [5] After research, rename the clusters according to immune cell type : #####

TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'Mono',
                                   '1' = 'CD4 helper',
                                   '2' = 'Neutrophils',
                                   '3' = 'T cell',
                                   '4' = 'Mono',
                                   '5' = 'CD8 T',
                                   '6' = 'B',
                                   '7' = 'CD8 T',
                                   '8' = 'DCs',
                                   '9' = 'NKT',
                                   '10' = 'Tregs',
                                   '11' = 'NK',
                                   '12' = 'Platelets',
                                   '13' = 'Unknown',
                                   '14' = 'DCs', 
                                   '15' = 'Naive T',
                                   '16' = 'Naive T',
                                   '17' = 'pDCs')

saveRDS(treatment, "honours/results/FinalIndex/adjtreatment.rds")
saveRDS(TreatmentAnnotated, "honours/results/FinalIndex/TreatmentAnnotatedforLabels.rds")
TreatmentAnnotated <- read_rds("honours/results/FinalIndex/TreatmentAnnotatedforLabels.rds")
treatment <- read_rds("honours/results/FinalIndex/adjtreatment.rds")

##### [6] Visualisation with annotations #####
# UMAP plot with labels 
Idents(TreatmentAnnotated)
# Mono CD4_helper neutrophils T_cell CD8_T B DCs Tregs NK platelets unknown Naive_T pDCs


TreatmentAnnotated@meta.data$cell_type <- TreatmentAnnotated@active.ident

palette.b <- c("#15c284", #0 mono
             "#d72554", #1 CD4_helper
             "#7ac745", #2 neutrophils
             "#7e549f", #3 T_cell
             "#9b78c2", #4 CD8_T
            #"#fb836f", #5 CD8_T
             "#a0d9e9", #6 
             #"#9b78c2", #7 CD8_T
             "#6ab5ba", #8 DCs
             "#93aff5", #9 NK
             "#c674bc", #10 Tregs
           # "#81cfff", #11 NK
             "#8edecf", #12 platelets
             "#69a923", #13 unknown
            #  "#00a68e", #14 DCs
             "#d73f3f", #15 naive T 
           # "white", #16 naive T 
             "#a0d9e9") #17 pDC


p2 <- DimPlot(TreatmentAnnotated, 
              reduction = "umap", 
              pt.size = 1.5,
              label = TRUE, 
              label.color = "white",
              label.box = TRUE,
              label.size = 5, 
              repel = TRUE,
              group.by = "celltype") +  # Map to the cluster variable
  scale_color_manual(values = palette.b) +
  scale_fill_manual(values = palette.b) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "", title = "") 

?DimPlot()


# view conserved cell markers across conditions : 
# Visualize gene expression changes using a variation of feature plots ; 

FeaturePlot(treatment, 
            features = c("FCGR3A"), 
            split.by = "stim", 
            max.cutoff = 3, 
            cols = c("grey","red"))


# See what genes change in different conditions for cells of the same type : 
# The code is specifically looking for genes that change in expression consistently across all cell types of the same type (need to specify to look at particular cell type)

treatment$celltype.stim <- paste(Idents(treatment), treatment$stim, sep = "_")
treatment$celltype <- Idents(treatment)
Idents(treatment) <- "celltype.stim"
response <- FindMarkers(treatment, ident.1 = "A_STIM", ident.2 = "CTRL", verbose = FALSE)
response <- write.xlsx(response, "honours/results/DifExpGenes")


# Identify differentially expressed gene across conditions 
# comparative analyses between stimulated & control cells  : plot average expression in stim & control cells and look for outliers 

theme_set(theme_cowplot())

# 
CellType1 <- subset(treatment, 
                    idents = "CellType1")
Idents(CellType1) <- "stim"
ExpressionCellType1 <- as.data.frame(log1p(AverageExpression(CellType1, verbose = FALSE)$RNA))
ExpressionCellType1$gene <- rownames(ExpressionCellType1)

# plots of average gene expression of all genes in a cell type from controlled condition to stimulated condition 
dp1 <- gglot(ExpressionCellType1, 
             aes(CTRL, STIM)) + 
  geom() + 
  ggtitle("Cell Type 1")

# labelling the outliers (these are point where expression of a gene changed between control & stim)
dp1 <- LabelPoints(plot = dp1, 
                   #point = , 
                   repel = TRUE)


