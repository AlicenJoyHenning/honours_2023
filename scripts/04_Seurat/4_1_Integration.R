# INTEGRATION AND DIMENSIONALITY REDUCTION 
# Seurat workflow after preprocessing

##### [1] Load dependencies & datasets #####
getwd()

library(BiocManager)
BiocManager::install('limma') # BiocManager::install('limma')
BiocManager::install("Seurat", version = "4.1")
BiocManager::install('multtest')
BiocManager::install('metap')
BiocManager::install('xlsx')
BiocManager::install('pheatmap')
BiocManager::install('XLConnect')
BiocManager::install('Matrix')
BiocManager::install('Seurat')
install.packages('Matrix')
install.packages('installr')
install.Rtools()
install.packages('Matrix')
install.packages('Seurat')
remotes::install_version("Seurat", version = "4.1.1")

library(remotes)
library(installr)
library(dplyr)
library(writexl)
library(openxlsx)
library(XLConnect)
library(SeuratObject)
library(ggplot2)
library(grid)
library(Seurat)
library(tidyverse)
library(celldex)
library(patchwork)
library(readr)
library(Matrix)
library(metap)
library(openxlsx)
library(readxl)
library(xlsx)

remotes::install_version("Matrix", version = "1.5.3")

packageVersion("Matrix")
remove.packages("SeuratObject")

# Load datasets: alpha, lambda, and untreated using ReadMtx function : 
Amatrix <- "honours/work/1109/alpha/matrix.mtx.gz"
Abarcodes <- "honours/work/1109/alpha/barcodes.tsv.gz"
Afeatures <- "honours/work/1109/alpha/features.tsv.gz"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=0)
saveRDS(alpha, "honours/work/1109/alpha/alphaCountMatrix.rds")

Lmatrix <- "honours/work/1109/lambda/matrix.mtx.gz"
Lbarcodes <- "honours/work/1109/lambda/barcodes.tsv.gz"
Lfeatures <- "honours/work/1109/lambda/features.tsv.gz"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features = 0)
saveRDS(lambda, "honours/work/1109/lambda/lambdaCountMatrix.rds")

Umatrix <- "honours/work/1109/untreated/matrix.mtx.gz"
Ubarcodes <- "honours/work/1109/untreated/barcodes.tsv.gz"
Ufeatures <- "honours/work/1109/untreated/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features = 0)
saveRDS(untreated, "honours/work/1109/untreated/untreatedCountMatrix.rds")

# noted : sizes alpha, lambda, untreated : 193, 219, 243 MB
# fixes new sizes : 193,  219, 198 MB
# final sizes : 249 286 257 MB

# Alternative way to read in datasets using Read10X function : 
alpha <- Read10X(data.dir = "honours/work/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha, project='ifnalpha', min.cells=3, min.features=200)

# Quick access to save the worked on Seurat objects : 
alpha <- saveRDS(alpha, "honours/work/RObjects/")
lambda <- saveRDS(lambda, "honours/work/RObjects/")
untreated <- saveRDS(untreated, "honours/work/RObjects/")
treatment <- saveRDS(treatment, "honours/work/1109/treatment.rds")

# To read in the saved Seurat objects : 
treatment <- readRDS("honours/work/1109/treatment.rds")


##### [2] Perform cell and gene level quality control independently on the datasets #####

# Removing unwanted cells based on # genes expressed and mitochondrial gene expression
# Create meta data column for mitochondrial gene identified by genes starting with MT : 
alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

alpha[["percent.xist"]] <- PercentageFeatureSet(alpha, pattern = "XIST")
lambda[["percent.xist"]] <- PercentageFeatureSet(lambda, pattern = "XIST")
untreated[["percent.xist"]] <- PercentageFeatureSet(untreated, pattern = "XIST")


# Subset the seurat objects based on mitochondrial percentage and the number of features : 
alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

# Normalise according to LOG scale (good for integration) : 
# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
lambda <- NormalizeData(lambda, normalization.method = "LogNormalize", scale.factor = 10000)
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)
# 2 : scaling the data :
# The results of this are stored in object[["RNA"]]@scale.data
all.alpha.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.alpha.genes)
all.lambda.genes <- rownames(lambda)
lambda <- ScaleData(lambda, features = all.lambda.genes)
all.untreated.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.untreated.genes)

# 3 : feature selection (check dim(object@assays$RNA@counts)[2] to view genes and length(object@assays$RNA@var.features) to view 2000 variable genes)
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
lambda <- FindVariableFeatures(lambda, selection.method = "vst", nfeatures = 2000)
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)

saveRDS(alpha, "honours/work/1109/lambda/alphaCountMatrixFiltered.rds")
saveRDS(lambda, "honours/work/1109/lambda/lambdaCountMatrixFiltered.rds")
saveRDS(untreated, "honours/work/1109/lambda/untreatedCountMatrixFiltered.rds")


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
?IntegrateData()

##### [4] Perform integration #####

# Identify integration anchors : 

anchors <- FindIntegrationAnchors(
  object.list = TreatmentList, 
  reference = 3, # untreated
  anchor.features = TreatmentFeatures,
  normalization.method = "LogNormalize"
)
# Darisia Index : Found 13669 anchors, retained 
# Index with problems; Found 12273 anchors, retained 5914
# Final Index : Found anchors 8220, retained 4846


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
colnames(treatment@meta.data)[1] <- "treatment"

# View whether the treatment datasets cluster together : 
colours <- c("#c35cad","#6ab5ba","#d3d3d3")
p1 <- DimPlot(treatment, reduction = "umap",pt.size = 1, group.by = "treatment") + scale_color_manual(values = colours) + ggtitle("")

# View the total number of clusters :
# Define color palette : 
palette <- c("#15c284", #0
               "#d72554", #1
               "#7ac745", #2
               "#7e549f", #3
               "#37a777", #4
               "#fb836f", #5 
               "#a0d9e9", #6
               "#9b78c2", #7
               "#6ab5ba", #8
               "#93aff5", #9
               "#c674bc", #10
               "#81cfff", #11
               "#8edecf", #12
               "#69a923", #13 
               "#00a68e", #14 
               "#d73f3f", #15
               "white", #16
               "#a0d9e9" #17
)
p2 <- DimPlot(treatment, reduction = "umap", pt.size = 1.5, label = TRUE, label.color = "black",label.box = FALSE, label.size = 6, repel = TRUE) + scale_color_manual(values = palette)
