# INTEGRATION AND DIMENSIONALITY REDUCTION 
# Seurat v3 workflow for pre-processing, dimensionality reduction, and visualisation 

##### [1] Load dependencies #####

#  Package names 
packages <- c("BiocVersion", "Seurat","Matrix","writexl","openxlsx","readxl","dplyr","tidyverse")

# Install packages not yet installed 
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

#  Packages loading
invisible(lapply(packages, library, character.only = TRUE))


##### [2] Load datasets using Read10X function ####

# Access folders in directory where matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz are stored
# Windows home computer : honours/work/1109/alpha/" 
# Ubuntu lab compute : /home/alicen/2024/interferon_dataset/alpha/"

getwd()

alpha <- Read10X("/home/alicen/2024/interferon_dataset/alpha/")
alpha <- CreateSeuratObject(alpha, project="alpha",  min.cells=3, min.features=0)

lambda <- Read10X("/home/alicen/2024/interferon_dataset/lambda/")
lambda <- CreateSeuratObject(lambda, project="lambda", min.cells=3, min.features = 0)

untreated <- Read10X("/home/alicen/2024/interferon_dataset/untreated/")
untreated <- CreateSeuratObject(untreated, project="untreated", min.cells=3, min.features = 0)


##### [3] Perform cell and gene level quality control independently on the datasets #####

# Removing unwanted cells based on # genes expressed and mitochondrial gene expression

# Create meta data column for mitochondrial gene identified by genes starting with MT : 
alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

# Subset the Seurat objects based on mitochondrial percentage and the number of features : 
alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)


# Normalize according to LOG scale (good for integration) : 
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


saveRDS(alpha, "/home/alicen/2024/interferon_dataset/R_objects/alphaCountMatrixFiltered.rds")
saveRDS(lambda, "/home/alicen/2024/interferon_dataset/R_objects/lambdaCountMatrixFiltered.rds")
saveRDS(untreated, "/home/alicen/2024/interferon_dataset/R_objects/untreatedCountMatrixFiltered.rds")


##### [4] Perform integration among the datasets #####

# Create a list of Seurat objects : 
TreatmentList <- list(alpha, lambda, untreated) 
# Select features that are repeatedly variable across datasets for integration : 
TreatmentFeatures <- SelectIntegrationFeatures(object.list = TreatmentList)

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
DefaultAssay(treatment) <- "integrated"  # designating an integrated data assay for a Seurat object. This step is crucial for integrating multiple scRNA-seq datasets and performing downstream analyses on the integrated data

##### [5] Dimensionality reduction #####

# Run standard workflow for visualization and clustering : 
treatment <- ScaleData(treatment, verbose = FALSE)
treatment <- RunPCA(treatment, npcs = 30, verbose = FALSE)
treatment <- RunUMAP(treatment, reduction = "pca", dims = 1:30) # 30
treatment <- FindNeighbors(treatment, reduction = "pca", dims = 1:30) # 30 
treatment <- FindClusters(treatment, resolution = 0.5)


##### [6] Preliminary visualization of clusters #####
# Change name of metadata column header to stim : 
colnames(treatment@meta.data)[1] <- "treatment"

# View whether the treatment datasets cluster together : 
colours <- c("#c35cad","#6ab5ba","#d3d3d3")
p1 <- DimPlot(treatment, reduction = "umap",pt.size = 1,group.by = "treatment") + scale_color_manual(values = colours) + ggtitle("") 
p1

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
p2



# Quick access to save the worked on Seurat objects ####
alpha <- saveRDS(alpha, "honours/work/RObjects/")
lambda <- saveRDS(lambda, "honours/work/RObjects/")
untreated <- saveRDS(untreated, "honours/work/RObjects/")
treatment <- saveRDS(treatment, "honours/work/1109/treatment.rds")

# To read in the saved Seurat objects : 
treatment <- readRDS("honours/work/1109/treatment.rds")

# Working on Loop to preprocess individual datasets #####
# Loop through the Seurat objects instead 
projects <- list(alpha, lambda, untreated)

for (project in projects) {
  print("Processing project:", project, "\n")
  # Create meta data column for mitochondrial genes
  project[["percent.mt"]] <- PercentageFeatureSet(get(project), pattern = "^MT-")
  project <- subset(project, project[["nFeature_RNA"]] > 200 & project[["nFeature_RNA"]] < 2500 & percent.mt < 10) # have to reference the metadata column specifically or else it won't be noticed
  project <- subset(project,project[["nCount_RNA"]] > 10)
  print("Quality control complete","\n")
  # Normalize the data
  test <- class(project)
  print(test)
  project <- NormalizeData(project, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  print("Data normalized for project:", project, "\n")
}




Scale the data
all_genes <- rownames(project)
ScaleData(project, features = all_genes)
cat("Data scaled for project:", project, "\n")

# Find variable features
FindVariableFeatures(project, selection.method = "vst", nfeatures = 2000) 
cat("Variable features identified for project:", project, "\n")

# Save the filtered Seurat object
# saveRDS(project), paste("honours/work/1109/", project, "/", project, "CountMatrixFiltered.rds", sep = ""))

cat("Filtered Seurat object saved for project:", project, "\n\n")


+



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



# WORKING EDIT: change subsetting to be based on outliers not arbitrary boundaries
# (1) extract vector of all Feature values 
x <- alpha[["nFeature_RNA"]]

