# R Script : Seurat workflow for integration of datasets after preprocessing with SASCRiP (before normalization with sctransform)

##### Loading dependencies and datasets #####

library(Seurat)
library(SeuratData)
library(patchwork)

# Load datasets: alpha, lambda, and untreated
# can i do readRDS instead ?
alpha <- Read10X(data.dir = "path/to/alpha/data")  
lambda <- Read10X(data.dir = "path/to/lambda/data")
untreated <- Read10X(data.dir = "path/to/untreated/data")

# Create a list of Seurat objects
treatment.list <- list(alpha, lambda, untreated)

##### Normalization before integration #####

# Normalize, find variable features, and select integration features for each dataset

treatment.list <- lapply(
  datasets, # list containing 3 datasets 
  function(x) {
    x <- NormalizeData(x)  # normalizes each dataset independently 
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) # identifies variable features for each dataset independently 
    return(x)
}) # here our output is the same list containing the same seurat objects that themselves have been modified. 

##### Prepare for integration #####

# Select features that are repeatedly variable across datasets for integration : 
treatment.features <- SelectIntegrationFeatures(object.list = treatment.list)
# function does not directly modify the Seurat objects themselves. Instead, it processes the gene expression data within the Seurat objects to identify a set of integration features, which are then used for subsequent integration steps. The Seurat objects remain unchanged, but the integration features selected based on these objects play a crucial role in guiding the integration process.









##### perform integration #####

# Identify integration anchors 
anchors <- FindIntegrationAnchors(
  object.list = dataset.list, 
  anchor.features = features
)

# Integrate datasets
combined <- IntegrateData(anchorset = anchors)

# Specify integrated data assay
DefaultAssay(combined) <- "integrated"

# Run standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
