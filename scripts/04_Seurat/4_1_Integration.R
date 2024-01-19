# INTEGRATION AND DIMENSIONALITY REDUCTION 
# Seurat workflow after preprocessing

##### [1] Load dependencies #####
getwd()
library(BiocManager)
BiocManager::install("Seurat", version = "4.1")
BiocManager::install('Matrix')
BiocManager::install('xlsx')
BiocManager::install('XLConnect')
BiocManager::install('writexl')
BiocManager::install('openxlsx')
BiocManager::install('readxl')
BiocManager::install('dplyr')

# data manipulation ####
library(tidyverse)
library(dplyr)
# storing data outputs in excel sheet ####
library(writexl) 
library(openxlsx)
library(XLConnect)
library(readxl)
# Seurat package (we love her) ####
library(SeuratObject) 
library(Seurat)
# plotting help ####
library(ggplot2) 
library(grid)
library(patchwork)
# support ####
library(readr)
library(Matrix)
library(metap)

# VIEWING PACKAGE VERSION #s for DEVICE COMPARISONS : 
# remotes::install_version("Matrix", version = "1.5.3")
# packageVersion("Matrix")
# remove.packages("SeuratObject")


# Load datasets: alpha, lambda, and untreated using Read10X function : ####
alpha <- Read10X("honours/work/1109/alpha/")
alpha <- CreateSeuratObject(alpha, project="alpha",  min.cells=3, min.features=0)

lambda <- Read10X("honours/work/1109/lambda/")
lambda <- CreateSeuratObject(lambda, project="lambda", min.cells=3, min.features = 0)

untreated <- Read10X("honours/work/1109/untreated/")
untreated <- CreateSeuratObject(untreated, project="untreated", min.cells=3, min.features = 0)





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

# WORKING EDIT: change subsetting to be based on outliers not arbitrary boundaries ####
# (1) extract vector of all Feature values 
x <- alpha("nFeature_RNA"]]



# Normalise according to LOG scale (good for integration) : ####
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

##### [4] Perform integration #####

# Identify integration anchors : 

anchors <- FindIntegrationAnchors(
  object.list = TreatmentList, 
  reference = 3, # untreated
  anchor.features = TreatmentFeatures,
  normalization.method = "LogNormalize"
)

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

##### Working on Loop to preprocess individual datasets #####
# Define the project names and data directories :
projects <- c("alpha", "lambda", "untreated")
data_dirs <- c("honours/work/1109/alpha", "honours/work/1109/lambda", "honours/work/1109/untreated")

# Loop through each dataset :
for (i in 1:length(projects)) {
  
  project_name <- projects[i]
  data_dir <- data_dirs[i]
  cat("Loading dataset for project:", project_name, "\n")
  
  # Read in the data using Read10X Seurat function & create Seurat object 
  project <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = project, project = project_name, min.cells = 0, min.features = 0)
  cat(project_name,"Seurat object created","\n")
  
  # Create meta data column for mitochondrial genes
  percent.mt <- paste("percent.mt", project_name, sep = ".")
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") 
  
  # Subset based on mitochondrial percentage and the number of features
  nFeature_RNA <- seurat_obj@meta.data$nFeature_RNA
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  cat("Filtering of", project_name,"complete.", "\n")
  
  # Log normalize the dataset
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  cat("Data normalized for project:", project_name, "\n")
  
  # Scale the data
  all_genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all_genes)
  cat("Data scaled for project:", project_name, "\n")
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  cat("Variable features identified for project:", project_name, "\n")
  
  # Save the filtered Seurat object
  saveRDS(seurat_obj, paste("honours/work/1109/", project_name, "/", project_name, "CountMatrixFiltered.rds", sep = ""))
  
  cat("", project_name, "\n\n")
}


# Define the file paths for each dataset
dataset_paths <- list(
  list(
    matrix = "honours/work/1109/alpha/matrix.mtx.gz",
    barcodes = "honours/work/1109/alpha/barcodes.tsv.gz",
    features = "honours/work/1109/alpha/features.tsv.gz",
    project = "alpha"
  ),
  list(
    matrix = "honours/work/1109/lambda/matrix.mtx.gz",
    barcodes = "honours/work/1109/lambda/barcodes.tsv.gz",
    features = "honours/work/1109/lambda/features.tsv.gz",
    project = "lambda"
  ),
  list(
    matrix = "honours/work/1109/untreated/matrix.mtx.gz",
    barcodes = "honours/work/1109/untreated/barcodes.tsv.gz",
    features = "honours/work/1109/untreated/features.tsv.gz",
    project = "untreated"
  )
)

# Loop through each dataset
for (dataset_info in dataset_paths) {
  cat("Loading dataset for project:", dataset_info$project, "\n")
  
  # Read in the matrix, barcodes, and features
  matrix_file <- dataset_info$matrix
  barcodes_file <- dataset_info$barcodes
  features_file <- dataset_info$features
  
  AlphaMatrix <- ReadMtx(matrix_file, barcodes_file, features_file)
  
  # Create a Seurat object
  project_name <- dataset_info$project
  assign(project_name, CreateSeuratObject(AlphaMatrix, project = project_name, min.cells = 3, min.features = 0))
  
  cat("Seurat object created for project:", dataset_info$project, "\n")
}

# Loop through the Seurat objects and perform the remaining steps
projects <- c("alpha", "lambda", "untreated")

for (project in projects) {
  cat("Processing project:", project, "\n")
  # Create meta data column for mitochondrial genes
  project[["percent.mt"]] <- PercentageFeatureSet(get(project), pattern = "^MT-")
}


cat("Mitochondrial gene percentages calculated for project:", project, "\n")

# Subset based on mitochondrial percentage and the number of features
subset(project, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent_mt_col < 10)
cat("Seurat object subset for project:", project, "\n")

# Normalize the data
NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
cat("Data normalized for project:", project, "\n")

# Scale the data
all_genes <- rownames(project)
ScaleData(project, features = all_genes)
cat("Data scaled for project:", project, "\n")

# Find variable features
FindVariableFeatures(project, selection.method = "vst", nfeatures = 2000) 
cat("Variable features identified for project:", project, "\n")

# Save the filtered Seurat object
# saveRDS(project), paste("honours/work/1109/", project, "/", project, "CountMatrixFiltered.rds", sep = ""))

cat("Filtered Seurat object saved for project:", project, "\n\n")
}





for (project in projects) {
  cat("Processing project:", project, "\n")
  
  # Create meta data column for mitochondrial genes
  percent_mt_col <- paste("percent.mt", project, sep = ".")
  assign(project, PercentageFeatureSet(get(project), pattern = "^MT-"))
  
  cat("Mitochondrial gene percentages calculated for project:", project, "\n")
  
  # Subset based on mitochondrial percentage and the number of features
  subset_expr <- paste("nFeature_RNA > 200 & nFeature_RNA < 2500 &", percent_mt_col, "< 10")
  assign(project, subset(get(project), subset = eval(parse(text = subset_expr))))
  cat("Seurat object subset for project:", project, "\n")
  
  # Normalize the data
  assign(project, NormalizeData(get(project), normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE))
  cat("Data normalized for project:", project, "\n")
  
  # Scale the data
  all_genes <- rownames(get(project))
  assign(project, ScaleData(get(project), features = all_genes))
  cat("Data scaled for project:", project, "\n")
  
  # Find variable features
  assign(project, FindVariableFeatures(get(project), selection.method = "vst", nfeatures = 2000)) 
  cat("Variable features identified for project:", project, "\n")
  
  # Save the filtered Seurat object
  # saveRDS(get(project), paste("honours/work/1109/", project, "/", project, "CountMatrixFiltered.rds", sep = ""))
  
  cat("Filtered Seurat object saved for project:", project, "\n\n")
}


##### 


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
