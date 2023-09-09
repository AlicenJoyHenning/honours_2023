# MANUAL CLUSTER ANNOTATION 
# Seurat workflow after integration

##### [1] Load dependencies $ datasets #####
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
library(writexl)
library(openxlsx)
library(XLConnect)
library(ggplot2)
library(grid)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(SingleR)
library(celldex)
library(patchwork)
library(readr)
library(Matrix)
library(metap)
library(openxlsx)
library(cowplot)
library(readxl)
library(xlsx)

# To read in the saved Seurat objects : 
treatment <- readRDS("honours/results/integrated.trials/treatmentsucess.rds")


##### [2] Identify Conserved cell type markers for cell type annotation #####

# Find markers for each cluster :
DefaultAssay(treatment) <- "RNA"

ClusterPath <- "honours/results/IntegratedMarkers/SeuratClusters.xlsx"

Cluster0Markers <- FindConservedMarkers(treatment, ident.1= 0, grouping.var = "stim", verbose = FALSE)
Cluster0Markers <- subset(Cluster0Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster0Markers <- subset(Cluster0Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C0Rows <- rownames(Cluster0Markers)
Cluster0Markers <- cbind(RowNames = C0Rows, Cluster0Markers)
write.xlsx(Cluster0Markers, ClusterPath, sheetName = "Cluster0")
dim(Cluster0Markers) # Before adjustments 1755 > 1217 > 57


Cluster1Markers <- FindConservedMarkers(treatment, ident.1= 1, grouping.var = "stim", verbose = FALSE)
Cluster1Markers <- subset(Cluster1Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster1Markers <- subset(Cluster1Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C1Rows <- rownames(Cluster1Markers)
Cluster1Markers <- cbind(RowNames = C1Rows, Cluster1Markers)
openxlsx::write.xlsx(Cluster1Markers, ClusterPath, sheetName = "Cluster1", append = TRUE)
dim(Cluster1Markers) # Before adjustments : 1841 > 1368 > 47

Cluster2Markers <- FindConservedMarkers(treatment, ident.1= 2, grouping.var = "stim", verbose = FALSE)
Cluster2Markers <- subset(Cluster2Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster2Markers <- subset(Cluster2Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C2Rows <- rownames(Cluster2Markers)
Cluster2Markers <- cbind(RowNames = C2Rows, Cluster2Markers)
dim(Cluster2Markers) # Before adjustments : 1730 > 99


Cluster2v3Markers <- FindConservedMarkers(treatment, ident.1= 2,ident.2 = c(3, 5, 7, 8, 9, 11, 14), grouping.var = "stim", verbose = FALSE)
Cluster2v3Markers <- subset(Cluster2v3Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster2v3Markers <- subset(Cluster2v3Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C2v3Rows <- rownames(Cluster2v3Markers)
Cluster2v3Markers <- cbind(RowNames = C2v3Rows, Cluster2v3Markers)
dim(Cluster2v3Markers) # 411 > 387 > 0


Cluster3Markers <- FindConservedMarkers(treatment, ident.1= 3, grouping.var = "stim", verbose = FALSE)
Cluster3Markers <- subset(Cluster3Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster3Markers <- subset(Cluster3Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C3Rows <- rownames(Cluster3Markers)
Cluster3Markers <- cbind(RowNames = C3Rows, Cluster3Markers)
dim(Cluster3Markers) # Before adjustments : 1637 > 1043 > 26

Cluster4Markers <- FindConservedMarkers(treatment, ident.1= 4, grouping.var = "stim", verbose = FALSE)
Cluster4Markers <- subset(Cluster4Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster4Markers <- subset(Cluster4Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C4Rows <- rownames(Cluster4Markers)
Cluster4Markers <- cbind(RowNames = C4Rows, Cluster4Markers)
dim(Cluster4Markers) # 1699  > 1237 > 46

Cluster5Markers <- FindConservedMarkers(treatment, ident.1= 5, grouping.var = "stim", verbose = FALSE)
Cluster5Markers <- subset(Cluster5Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster5Markers <- subset(Cluster5Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C5Rows <- rownames(Cluster5Markers)
Cluster5Markers <- cbind(RowNames = C5Rows, Cluster5Markers)
dim(Cluster5Markers) # 1525 > 821 > 57

Cluster6Markers <- FindConservedMarkers(treatment, ident.1= 6, grouping.var = "stim", verbose = FALSE)
Cluster6Markers <- subset(Cluster6Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster6Markers <- subset(Cluster6Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C6Rows <- rownames(Cluster6Markers)
Cluster6Markers <- cbind(RowNames = C6Rows, Cluster6Markers)
dim(Cluster6Markers) # 1743 > 1023 > 52


Cluster7Markers <- FindConservedMarkers(treatment, ident.1= 7, grouping.var = "stim", verbose = FALSE)
Cluster7Markers <- subset(Cluster7Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster7Markers <- subset(Cluster7Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C7Rows <- rownames(Cluster7Markers)
Cluster7Markers <- cbind(RowNames = C7Rows, Cluster7Markers)
dim(Cluster7Markers) # 1576 > 845 > 39

Cluster8Markers <- FindConservedMarkers(treatment, ident.1= 8, grouping.var = "stim", verbose = FALSE)
Cluster8Markers <- subset(Cluster8Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster8Markers <- subset(Cluster8Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C8Rows <- rownames(Cluster8Markers)
Cluster8Markers <- cbind(RowNames = C8Rows, Cluster8Markers)
dim(Cluster8Markers) # 1592 > 871 > 37

Cluster9Markers <- FindConservedMarkers(treatment, ident.1= 9, grouping.var = "stim", verbose = FALSE)
Cluster9Markers <- subset(Cluster9Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster9Markers <- subset(Cluster9Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C9Rows <- rownames(Cluster9Markers)
Cluster9Markers <- cbind(RowNames = C9Rows, Cluster9Markers)
dim(Cluster9Markers) # 1743 > 740 > 77

Cluster10Markers <- FindConservedMarkers(treatment, ident.1= 10, grouping.var = "stim", verbose = FALSE)
Cluster10Markers <- subset(Cluster10Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster10Markers <- subset(Cluster10Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C10Rows <- rownames(Cluster10Markers)
Cluster10Markers <- cbind(RowNames = C10Rows, Cluster10Markers)
dim(Cluster10Markers) # 1247 > 430 > 100 

Cluster11Markers <- FindConservedMarkers(treatment, ident.1= 11, grouping.var = "stim", verbose = FALSE)
Cluster11Markers <- subset(Cluster11Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster11Markers <- subset(Cluster11Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C11Rows <- rownames(Cluster11Markers)
Cluster11Markers <- cbind(RowNames = C11Rows, Cluster11Markers)
dim(Cluster11Markers) #1388 > 384 > 10 

Cluster12Markers <- FindConservedMarkers(treatment, ident.1= 12, grouping.var = "stim", verbose = FALSE)
Cluster12Markers <- subset(Cluster12Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster12Markers <- subset(Cluster12Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C12Rows <- rownames(Cluster12Markers)
Cluster12Markers <- cbind(RowNames = C12Rows, Cluster12Markers)
dim(Cluster12Markers) # 1847 > 199 > 127 

Cluster13Markers <- FindConservedMarkers(treatment, ident.1= 13, grouping.var = "stim", verbose = FALSE)
Cluster13Markers <- subset(Cluster13Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster13Markers <- subset(Cluster13Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C13Rows <- rownames(Cluster13Markers)
Cluster13Markers <- cbind(RowNames = C13Rows, Cluster13Markers)
dim(Cluster13Markers) # 1567 > 58 > 54 

Cluster14Markers <- FindConservedMarkers(treatment, ident.1= 14, grouping.var = "stim", verbose = FALSE)
Cluster14Markers <- subset(Cluster14Markers, subset = lambda_p_val_adj < 0.05 & alpha_p_val_adj < 0.05 & untreated_p_val_adj < 0.05) 
Cluster14Markers <- subset(Cluster14Markers, subset = lambda_avg_log2FC > 1 & alpha_avg_log2FC > 1 & untreated_avg_log2FC > 1) 
C14Rows <- rownames(Cluster14Markers)
Cluster14Markers <- cbind(RowNames = C14Rows, Cluster14Markers)
dim(Cluster14Markers) # 1355 > 180 > 79 


# Saving the data to an excel sheet to annotate genes : 
clusters <- list('Cluster0' = Cluster0Markers, 'Cluster1' = Cluster1Markers,'Cluster2' = Cluster2Markers, 'Cluster3' = Cluster3Markers,'Cluster4' = Cluster4Markers, 'Cluster5' = Cluster5Markers, 'Cluster6' = Cluster6Markers, 'Cluster7' = Cluster7Markers, 'Cluster8' = Cluster8Markers, 'Cluster9' = Cluster9Markers, 'Cluster10' = Cluster10Markers, 'Cluster11' = Cluster11Markers,'Cluster12' = Cluster12Markers,'Cluster13' = Cluster13Markers,'Cluster14' = Cluster14Markers)
openxlsx::write.xlsx(clusters, file = "honours/results/IntegratedMarkers/SeuratMarkers.xlsx")
openxlsx::write.xlsx(Cluster2v3Markers, file = "honours/results/IntegratedMarkers/Cluster2v3Markers.xlsx")

##### [3] Tools for Manual Annotation checking literature markers  #####

#Cluster #
features = c("CD41", "CD42b", "CD62P", "CD63")

FeaturePlot(object = treatment, 
            features = features,
            cols = c("lightgrey", "black"),
            label = TRUE,
            pt.size = 1.5, 
            blend = FALSE, 
            interactive = FALSE) + theme(
              panel.background = element_rect(fill = "darkgrey")) 

cluster_boundaries <- c(5, 7)

DotPlot(object = treatment, 
        features = features,
        cols = c("grey", "#265221")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Adjust x-axis labels angle
  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

##### [4] After research, rename the clusters according to immune cell type : #####

TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'Mono',
                                   '1' = 'Mono',
                                   '2' = 'CD4_helper',
                                   '3' = 'CD4_naive',
                                   '4' = 'Neutro',
                                   '5' = 'CD8',
                                   '6' = 'B',
                                   '7' = 'CD8',
                                   '8' = 'CD8',
                                   '9' = 'NK',
                                   '10' = '10',
                                   '11' = 'Tregs',
                                   '12' = '12',
                                   '13' = '13',
                                   '14' = '14')

saveRDS(treatment, "honours/results/IntegratedMarkers/treatment.rds")
saveRDS(TreatmentAnnotated, "honours/results/IntegratedMarkers/TreatmentAnnotated.rds")

##### [5] Visualisation with annotations #####
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


