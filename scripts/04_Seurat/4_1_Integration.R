# INTEGRATION AND DIMENSIONALITY REDUCTION 
# Seurat workflow after preprocessing

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

# Load datasets: alpha, lambda, and untreated using ReadMtx function : 
AlphaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/matrix.mtx.gz","honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)

LambdaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)

UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz")
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)



# noted : sizes alpha, lambda, untreated : 193, 219, 243 MB
# fixes new sizes : 193,  219, 198

# Alternative way to read in datasets using Read10X function : 
alpha <- Read10X(data.dir = "honours/work/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha, project='ifnalpha', min.cells=3, min.features=200)

# Quick access to save the worked on Eeurat objects : 
alpha <- saveRDS(alpha, "honours/work/RObjects/")
lambda <- saveRDS(lambda, "honours/work/RObjects/")
untreated <- saveRDS(untreated, "honours/work/RObjects/")
treatment <- saveRDS(treatment, "honours/results/integrated.trials/treatmentsucess.rds")

# To read in the saved Seurat objects : 
treatment <- readRDS("honours/results/integrated.trials/treatmentsucess.rds")

##### [2] Perform cell and gene level quality control independently on the datasets #####

# Removing unwanted cells based on # genes expressed and mitochondrial gene expression
# Create meta data column for mitochondrial gene identified by genes starting with MT : 
alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

# Subset the seurat objects based on mitochondrial percentage and the number of features : 
alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Normalise according to LOG scale (good for integration) : 
# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
lambda <- NormalizeData(lambda, normalization.method = "LogNormalize", scale.factor = 10000)
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)

# 2 : feature selection (check dim(object@assays$RNA@counts)[2] to view genes and length(oject@assays$RNA@var.features) to view 2000 variable genes)
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
lambda <- FindVariableFeatures(lambda, selection.method = "vst", nfeatures = 2000)
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)
# 
# # 3 : scaling the data :
# The results of this are stored in object[["RNA"]]@scale.data
all.alpha.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.alpha.genes)
all.lambda.genes <- rownames(lambda)
lambda <- ScaleData(lambda, features = all.lambda.genes)
all.untreated.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.untreated.genes)

##### [3] Prepare datasets for integration #####

# Create a list of Seurat objects : 
TreatmentList <- list(alpha, lambda, untreated) 

# Select features that are repeatedly variable across datasets for integration : 
TreatmentFeatures <- SelectIntegrationFeatures(object.list = TreatmentList) 
# Selecting features (genes) that are consistently variable across multiple
# datasets for the purpose of integrating those datasets into a single analysis
# function does not directly modify the Seurat objects instead, 
# processes the gene expression data within the Seurat objects to identify a set
# of integration features, which are then used for subsequent integration steps.
# The Seurat objects remain unchanged, but the integration features selected
# based on these objects play a crucial role in guiding the integration
# process.it goes as separate input into the integration function

?SelectIntegrationFeatures()
?FindIntegrationAnchors()

##### [4] Perform integration #####

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


# Run standard workflow for visualization and clustering : 
treatment <- ScaleData(treatment, verbose = FALSE)
treatment <- RunPCA(treatment, npcs = 30, verbose = FALSE)
treatment <- RunUMAP(treatment, reduction = "pca", dims = 1:30) # 30
treatment <- FindNeighbors(treatment, reduction = "pca", dims = 1:30) # 30 
treatment <- FindClusters(treatment, resolution = 0.5)

##### [5] Preliminary visualization of clusters #####
# Change name of metadata column header to stim : 
colnames(treatment@meta.data)[1] <- "stim"

# View whether the treatment datasets cluster together : 
colours <- c("#c35cad","#6ab5ba","#d3d3d3")
p1 <- DimPlot(treatment, reduction = "umap",pt.size = 1.5, group.by = "stim") + scale_color_manual(values = colours)

# View the total number of clusters :
# Define color palette : 
palette.b <- c("#FB836F", #0
               "#d72554", #1
               "#6ab5ba", #2
               "#2e8f95", #3
               "#7E549F", #4
               "#8caf2e", #5
               "#69a923", #6
               "#297b57", #7 
               "#00945a", #8
               "#265221", #9
               "#FFCB3E", #10
               "#00a68e", #11
               "#5c040c", #12
               "#ef931b", #13
               "#6ab5ba" #14
)
p2 <- DimPlot(treatment, reduction = "umap", pt.size = 1.5, label = TRUE, label.size = 6, repel = TRUE) + scale_color_manual(values = palette.b)

