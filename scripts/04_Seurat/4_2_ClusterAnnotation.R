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
treatment <- readRDS("honours/results/FinalIndex/adjtreatment.rds")

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
DefaultAssay(treatment)

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

# FeaturePlots 

FeaturePlot(treatment, features, cols = c("grey", "black"))

##### [4] Tools for Manual Annotation checking literature markers  #####

#Cluster #
DefaultAssay(treatment) <- "RNA"

features = c("SRY", "DDX3Y", "ZFY", "TSPY", "PRY", "UTY", "XIST")
features = c("CD8A", "ZNF683", "CD8B", "GZMB", "GZMA", "CD3G", "CCR7","CX3CR1", "DSTN")
features = c("CD3G", "IL17R", "CCR4", "CX3CR1", "CD161", "CCR6", "Leu8")
features = c("FCGR3A", "CD14") 
features <- ("CLEC4A")  #"CD4","CCR6", "KLF2", "IL6ST", 
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
# 
# TreatmentAnnotated <- RenameIdents(TreatmentAnnotated, 
#                                    'CD4+_helper' = 'CD4+ helper')

saveRDS(treatment, "honours/results/FinalIndex/adjtreatment.rds")
saveRDS(TreatmentAnnotated, "honours/results/FinalIndex/TreatmentAnnotated.rds")


##### [6] Visualisation with annotations #####

TreatmentAnnotated <- read_rds("honours/results/FinalIndex/TreatmentAnnotated.rds")

# UMAP plot with labels 
Idents(TreatmentAnnotated)
# Mono CD4_helper neutrophils T_cell CD8_T B DCs Tregs NK platelets unknown Naive_T pDCs


TreatmentAnnotated@meta.data$cell_type <- TreatmentAnnotated@active.ident
TreatmentAnnotated <- subset(TreatmentAnnotated, subset = seurat_clusters != 13) 

palette.b <- c("#15c284", #0 mono
               "#d72554", #1 naive CD4
               "#7ac745", #2 neutrophils
               "#9b78c2", #3 T helper
               # "#9b78c2", #4 mono
               "#dd5839", #5 naive CD8 T
               "#50b9d7", #6 B 
               "#936bb1", #7 cytotoxic T cells
               "#6ab5ba", #8 mDCs
               "#7e549f", #9 NKT
               "#c674bc", #10 Tregs
               "#709bfc", #11 NK
               "#47c9b0", #12 platelets
               "#92b832", #13 unknown
               "#00a68e", #14 DCs
               "#f46897", #15 & 16 memory T
               # "white", 
               "#55bbaa") #17 pDC


p2 <- DimPlot(TreatmentAnnotated, 
              reduction = "umap", 
              pt.size = 1.5,
              label = TRUE, 
              label.color = "white",
              label.box = TRUE,
              label.size = 7, 
              repel = TRUE,
              group.by = "cell_type") +  # Map to the cluster variable
  scale_color_manual(values = palette.b) +
  scale_fill_manual(values = palette.b) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2", color = "", title = "PBMC populations") 


features <- { c("CCL4L2",  # 0 & 4 : monocytes CCL4L2 IL1RN MNDA HCK
              "LEF1",  # 1 CD4 naive
              "CXCR2", # 2 neutrophils CXCR2
              "CCR4", # 3 T helper cells
              "CD8B",  # 5 naive CD8 T cells
              "MS4A1", # 6 B cells 
              "EOMES",  # 7 cytotoxic T  CCL5 EOMES
              "MAFB",  # 8 mDCs
              "KLRB1", # 9 NKT 
              "FOXP3", # 10 Treg
              "KLRF1", # 11 NK
              "GP9",   # 12 platelets (NRGN)
              "RFLNB", # 13 unknown 
              "VEGFA", # 14 DCs
              "KLF2",  # 15 & 16 memory T CCR7 KLF2 NRGN
              "MZB1")  # 17 pDcs 
}



overall <- DotPlot(object = TreatmentAnnotated, 
        features = "MS4A1",
        cols = c("white", "darkgrey")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
  scale_y_discrete(position = "left") +
 geom_abline(intercept = 0, slope = 1, color = "black",linetype = "dashed", size = 0.25)


BPlot <- DotPlot(object = TreatmentAnnotated, 
                 features = c("MS4A1", "MAFB", "MZB1", "GP9"),
                 cols = c("white", "#a9a9a9")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
  scale_y_discrete(position = "left") +
  geom_vline(xintercept = as.numeric(1.5, 3.5),  color = "black",linetype = "dashed", size = 0.25)

print(BPlot)


# T Cells 
TGroups <- c("NK", "NKT", "cytotoxic T", "naive CD8 T", "naive CD4 T", "T helper", "Tcm", "Tregs")
Tcells <- subset(TreatmentAnnotated, subset = cell_type %in% TGroups)
levels(Tcells)<- c("NK", "NKT", "cytotoxic T", "naive CD8 T", "Tregs", "naive CD4 T", "T helper", "Tcm")


TPlot <- DotPlot(object = Tcells, 
            features = c("CD3G","NKG7", "CD8A", "CD8B","FOXP3", "CD4","CD28", "LEF1"),
            cols = c("grey", "#a9a9a9")) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
      scale_y_discrete(position = "left") +
      geom_hline(yintercept = as.numeric(4.5),  color = "black",linetype = "dashed", size = 0.25)

print(TPlot)

# myeloid cells 
TreatmentM <- RenameIdents(treatment, 
                                   '0' = 'monocytes1',
                                   '2' = 'neutrophils',
                                   '4' = 'monocytes2',
                                   '8' = 'mDCs'
                                   )
TreatmentM@meta.data$cell_type <- TreatmentM@active.ident
MGroups <- c('monocytes1','monocytes2', 'neutrophils','mDCs')
Mcells <- subset(TreatmentM, subset = cell_type %in% MGroups)
levels(Mcells)<- MGroups


MPlot <- DotPlot(object = Mcells, 
                 features = c( "TLR6",  "CXCR1", "CXCR2", "FCGR3B", "ITGAX","TREM1","CD14", "FCGR3A"),
                 cols = c("white", "#a9a9a9")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
  scale_y_discrete(position = "left") +
  geom_hline(yintercept = as.numeric(4.5),  color = "black",linetype = "dashed", size = 0.25)

print(MPlot)

# Final total 
TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'monocytes(1)',
                                   '1' = 'naive CD4+ T',
                                   '2' = 'neutrophils',
                                   '3' = 'T helper',
                                   '4' = 'monocytes(2)',
                                   '5' = 'naive CD8+ T',
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

groups <- c('unknown', 'platelets', 'monocytes(1)', 'monocytes(2)', 'neutrophils', 
            'DCs', 'mDCs', 'pDCs', 
            'Tcm', 'T helper', 'Tregs', 'naive CD4+ T', 'naive CD8+ T', 'cytotoxic T', 'NKT', 'NK', 
            'B')
levels(TreatmentAnnotated) <- groups

DefaultAssay(TreatmentAnnotated) <- 'RNA'

Plot <- DotPlot(object = TreatmentAnnotated, 
                features = c("MS4A1", # B cells 
                             "CD3G","FOXP3","NKG7", "GNLY", "GZMB","KLRB1", "CD8A", "CD8B",  "LEF1", "CCR4", "GATA3", "CD4","SELL",  # T cells 
                  "TREM1","CD14","FCGR3B", "TLR6", "CXCR2", "SPI1",  # monocytes, neutrophils, 
                  "ITGAX","MAFB", "MZB1", # dendritic cells
                  "GP9" # platelets
                  ),
                cols = c("white", "darkgrey")) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.7)) +
  scale_y_discrete(position = "right") +
  geom_hline(yintercept = as.numeric(c(1.5, 2.5, 4.5, 5.5, 8.5, 16.5)),  color = "black", size = 0.2) +
  geom_vline(xintercept = as.numeric(c(1.5, 14.5, 20.5, 23.5)),  color = "black",linetype = "dashed", size = 0.2) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))

print(Plot)





##### [7] I can't remember what any of this is #####
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


