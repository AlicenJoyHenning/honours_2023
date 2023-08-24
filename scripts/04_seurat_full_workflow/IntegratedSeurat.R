# R Script : Seurat workflow for integration of datasets after preprocessing with SASCRiP (before normalization with sctransform)

##### Loading dependencies and datasets #####

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
#library(SeuratData)
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
library(wri)

# Load datasets: alpha, lambda, and untreated

# alpha <- Read10X(data.dir = "honours/work/ifnalpha/seurat_matrix/")
# alpha <- CreateSeuratObject(counts=alpha, project='ifnalpha', min.cells=3, min.features=200)
# lambda <- Read10X(data.dir = "honours/work/ifnlambda/seurat_matrix/")
# lambda <- CreateSeuratObject(counts=lambda, project='ifnlambda', min.cells=3, min.features=200)

# untreated object first requires the matrix to be made : 
UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz")
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)

AlphaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/matrix.mtx.gz","honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)

LambdaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)


# noted : sizes alpha, lambda, untreated : 193, 219, 243 MB
# fixes new sizes : 193,  219, 198

# matrix <- ReadMtx("honours/work/untreated/sm/seurat_matrix/matrix.mtx.gz", "honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz", "honours/work/untreated/sm/seurat_matrix/features.tsv.gz", skip.feature = 2)
# # Error: Matrix has 35639 rows but found 69537 features. Try increasing `skip.feature`. 

alpha <- saveRDS(alpha, "honours/work/RObjects/")
lambda <- saveRDS(lambda, "honours/work/RObjects/")
untreated <- saveRDS(untreated, "honours/work/RObjects/")
treatment <- saveRDS(treatment, "honours/results/integrated.trials/treatmentsucess.rds")


treatment <- readRDS("honours/results/integrated.trials/treatmentsucess.rds")

##### Perform quality control independently on the datasets #####

# Removing unwanted cells based on # genes expressed and mitochondrial gene expression

alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Normalising according to LOG scale (good for integration)
# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
lambda <- NormalizeData(lambda, normalization.method = "LogNormalize", scale.factor = 10000)
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)

# 2 : feature selection (check dim(object@assays$RNA@counts)[2] to view genes and length(oject@assays$RNA@var.features) to view 2000 variable genes)
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
lambda <- FindVariableFeatures(lambda, selection.method = "vst", nfeatures = 2000)
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)
# 
# # 3 : scaling the data : NO! actually doing this afterwards !!!
# # The results of this are stored in object[["RNA"]]@scale.data
# 
all.alpha.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.alpha.genes)
all.lambda.genes <- rownames(lambda)
lambda <- ScaleData(lambda, features = all.lambda.genes)
all.untreated.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.untreated.genes)

##### Prepare datasets for integration #####

TreatmentList <- list(alpha, lambda, untreated) # Create a list of Seurat objects

TreatmentFeatures <- SelectIntegrationFeatures(object.list = TreatmentList) # Select features that are repeatedly variable across datasets for integration : 
# electing features (genes) that are consistently variable across multiple datasets for the purpose of integrating those datasets into a single analysis
# function does not directly modify the Seurat objects themselves. Instead, it processes the gene expression data within the Seurat objects to identify a set of integration features, which are then used for subsequent integration steps. The Seurat objects remain unchanged, but the integration features selected based on these objects play a crucial role in guiding the integration process.
# it goes as separate input into the integration function

?SelectIntegrationFeatures()
?FindIntegrationAnchors()

##### perform integration #####

# Identify integration anchors : 

anchors <- FindIntegrationAnchors(
  object.list = TreatmentList, 
  reference = 3,
  anchor.features = TreatmentFeatures
)
# Found 13669 anchors, retained 


# Integrate data sets using the anchors : 
treatment <- IntegrateData(anchorset = anchors)

# Specify integrated data assay : 
DefaultAssay(treatment) <- "integrated"  
# designating an integrated data assay for a Seurat object. This step is crucial for integrating multiple scRNA-seq datasets and performing downstream analyses on the integrated data


# Run standard workflow for visualization and clustering
treatment <- ScaleData(treatment, verbose = FALSE)
treatment <- RunPCA(treatment, npcs = 30, verbose = FALSE)
treatment <- RunUMAP(treatment, reduction = "pca", dims = 1:30) # 30
treatment <- FindNeighbors(treatment, reduction = "pca", dims = 1:30) # 30 
treatment <- FindClusters(treatment, resolution = 0.5)

##### Visualization of clusters #####
# getting an error for 'stim' , clustered by seurat_clusters = wack, 
# changed name of metadata column header to stim
colnames(treatment@meta.data)[1] <- "stim"

library(viridisLite)

viridis <- viridisLite::viridis(15)

p1 <- DimPlot(treatment, reduction = "umap", group.by = "stim")
p2 <- DimPlot(treatment, reduction = "umap", pt.size = 1.5, label = TRUE, label.size = 6, label.box = TRUE, repel = TRUE)
# + scale_color_manual(values = viridis)
p3 <- DimPlot(treatment, reduction = "umap", split.by = "stim", label = TRUE)
p1 + p2



#### Identify Conserved cell type markers #####
# performing differential expression after integration :
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


###### Manual Annotation using literature markers  #####

#Cluster 0 
features = c("SOD2")
   
FeaturePlot(object = treatment, 
            features = features,
            cols = c("grey", "#225EA8"),
            label = TRUE,
            pt.size = 1.5, 
            blend = FALSE, 
            interactive = FALSE)
VlnPlot(treatment, features = "FCGR3A")
DotPlot()


##### Manual annotation confirmation using seurat markers #####

# check for validity of markers within this list of genes : 
# using the markers identified, a set of Feature plots will br made fo each gene (marker) to see if they acurately describe clusters: 


FeaturePlot <-  FeaturePlot(treatment, features = c("CCR7", "CD7"))
VlnPlot(treatment, features = "FCGR3A")
DotPlot()

###### After research, rename the clusters according to immune cell type : #####

TreatmentAnnotated <- RenameIdents(treatment, 
                          '0' = 'macrophages0',
                          '1' = 'macrophages1',
                          '2' = 'CD4+ T helper cells',
                          '3' = 'Naive CD4+ T cells',
                          '4' = 'monocytes',
                          '5' = 'CD8+ T cells',
                          '6' = 'B cells',
                          '7' = 'CD8+ T cells',
                          '8' = 'CD8+ T cells',
                          '9' = 'NK cells',
                          '10' = 'neutrophils',
                          '11' = 'Tregs',
                          '12' = 'DC',
                          '13' = 'macrophage?',
                          '14' = 'DC?')
saveRDS(treatment, "honours/results/IntegratedMarkers/treatment.rds")
saveRDS(TreatmentAnnotated, "honours/results/IntegratedMarkers/TreatmentAnnotated.rds")

p4 <- DimPlot(TreatmentAnnotated)
p5 <- DimPlot(TreatmentAnnotated, label = TRUE, label.box = TRUE)
p4+p5

##### VIsualisation with annotations #####
# view conserved cell markers across conditions : 
# expression level + pergentage of cells in a cluster 

Idents(treatment) <- factor(Idents(treatment), 
                            levels = c("identified cell types", "")
                            )
Markers <- c("gene markers to plot", "") # top 3 marker genes for each cluster 
features <- c("CD8B")
p5 <- DotPlot(treatment, 
              features,  
              cols = c("grey","#FF9AA2"), 
              dot.scale = 6)  #+ 
  #RotatedAxis()

p5 <- DotPlot(
  treatment,
  assay = NULL,
  features,
  cols = c("lightgrey", "#FF9AA2"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)


# See what genes change in different conditions for cells of the same type : 
# The code is specifically looking for genes that change in expression consistently across all cell types of the same type (need to specify to look at particular cell type)

treatment$celltype.stim <- paste(Idents(treatment), treatment$stim, sep = "_")
treatment$celltype <- Idents(treatment)
Idents(treatment) <- "celltype.stim"
response <- FindMarkers(treatment, ident.1 = "A_STIM", ident.2 = "CTRL", verbose = FALSE)
response <- write.xlsx(response, "honours/results/DifExpGenes")

# Visualize gene expression changes using a variation of feature plots ; 

FeaturePlot(treatment, 
            features = c("FCGR3A"), 
            split.by = "stim", 
            max.cutoff = 3, 
            cols = c("grey","red"))
# create fp showing gene expression in untreated vs treated cells (2 fp)


# Identify differentially expressed gene across conditions 
# comparative analyses btwn stim & control cells  : plot average expression in stim & control cells and look for outliers 

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



#### SingleR #####

reference <- BlueprintEncodeData()
reference <- ImmGenData() # didn't work 
reference <- DatabaseImmuneCellExpressionData()
reference <- MonacoImmuneData()
reference <- HumanPrimaryCellAtlasData()
celldex::MonacoImmuneData()

predictions <- SingleR(GetAssayData(treatment, assay = "RNA", slot = "data"),
                       clusters = Idents(treatment),
                       ref = reference,
                       labels = reference$label.fine)
Idents(treatment)
write.csv(predictions$labels,"honours/results/PredictionsBEDFine.txt")
plotScoreHeatmap(predictions)

SingleRBPETreatment <- RenameIdents(treatment,
                                    "0" = "Neutrophils",
                                    "1" = "CD4+ T-cells",
                                    "2" = "CD4+ T-cells",
                                    "3" = "Neutrophils",
                                    "4" = "Neutrophils",
                                    "5" = "CD4+ T-cells",
                                    "6" = "B-cells",
                                    "7" = "CD8+ T-cells",
                                    "8" = "NK cells",
                                    "9" = "CD4+ T-cells",
                                    "10" = "CD8+ T-cells",
                                    "11" = "Neutrophils",
                                    "12" = "Monocytes",
                                    "13" = "Neutrophils",
                                    "14" = "Eosinophils",
                                    "15" = "Monocytes")
                        
Idents(SingleRBPETreatment)
SingleRPlot <- DimPlot(SingleRBPETreatment,
                       reduction = "umap", 
                       label = TRUE)
SingleRPlot
