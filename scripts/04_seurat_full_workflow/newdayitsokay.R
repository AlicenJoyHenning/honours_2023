# R Script : Seurat workflow for integration of datasets after preprocessing with SASCRiP (before normalization with sctransform)

##### Loading dependencies and datasets #####

getwd()
library(BiocManager)
BiocManager::install("SeuratData")

library(dplyr)
library(ggplot2)
library(grid)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(patchwork)
library(readr)


# Load datasets: alpha, lambda, and untreated

alpha <- Read10X(data.dir = "honours/work/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha, project='ifnalpha', min.cells=3, min.features=200)
lambda <- Read10X(data.dir = "honours/work/ifnlambda/seurat_matrix/")
lambda <- CreateSeuratObject(counts=lambda, project='ifnlambda', min.cells=3, min.features=200)
untreated <- Read10X(data.dir = "honours/work/s/")
untreated <- CreateSeuratObject(counts=untreated, project='untreated', min.cells=3, min.features=200)
# noted : sizes alpha, lambda, untreated : 193, 219, 243 MB

# matrix <- ReadMtx("honours/work/untreated/sm/seurat_matrix/matrix.mtx.gz", "honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz", "honours/work/untreated/sm/seurat_matrix/features.tsv.gz", skip.feature = 2)
# # Error: Matrix has 35639 rows but found 69537 features. Try increasing `skip.feature`. 



alpha <- saveRDS(alpha, "honours/work/RObjects/")
lambda <- saveRDS(lambda, "honours/work/RObjects/")
untreated <- saveRDS(untreated, "honours/work/RObjects/")


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
# all.alpha.genes <- rownames(alpha)
# alpha <- ScaleData(alpha, features = all.alpha.genes)
# all.lambda.genes <- rownames(lambda)
# lambda <- ScaleData(lambda, features = all.lambda.genes)
# all.untreated.genes <- rownames(untreated)
# untreated <- ScaleData(untreated, features = all.untreated.genes)

##### Prepare datasets for integration #####

treatment.list <- list(alpha, lambda, untreated) # Create a list of Seurat objects

treatment.features <- SelectIntegrationFeatures(object.list = treatment.list) # Select features that are repeatedly variable across datasets for integration : 
# electing features (genes) that are consistently variable across multiple datasets for the purpose of integrating those datasets into a single analysis
# function does not directly modify the Seurat objects themselves. Instead, it processes the gene expression data within the Seurat objects to identify a set of integration features, which are then used for subsequent integration steps. The Seurat objects remain unchanged, but the integration features selected based on these objects play a crucial role in guiding the integration process.
# it goes as separate input into the integration function

##### perform integration #####

# Identify integration anchors : 

anchors <- FindIntegrationAnchors(
  object.list = treatment.list, 
  anchor.features = treatment.features
)

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

# Visualization
# getting an error for 'stim' , clustered by seurat_clusters = wack, 
# changed name of metadata column header to stim
colnames(treatment@meta.data)[1] <- "stim"


p1 <- DimPlot(treatment, reduction = "umap", group.by = "stim")
p2 <- DimPlot(treatment, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(treatment, reduction = "umap", split.by = "stim")
p1 + p2
