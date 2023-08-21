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

library(dplyr)
library()
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


# Load datasets: alpha, lambda, and untreated

# alpha <- Read10X(data.dir = "honours/work/ifnalpha/seurat_matrix/")
# alpha <- CreateSeuratObject(counts=alpha, project='ifnalpha', min.cells=3, min.features=200)
# lambda <- Read10X(data.dir = "honours/work/ifnlambda/seurat_matrix/")
# lambda <- CreateSeuratObject(counts=lambda, project='ifnlambda', min.cells=3, min.features=200)

# untreated object first requires the matrix to be made : 
UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", )
untreated <- CreateSeuratObject(matrix, project="untreated", min.cells=3, min.features=200)

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

TreatmentList <- list(alpha, lambda)# untreated) # Create a list of Seurat objects

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
  reference = 1,
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


p1 <- DimPlot(treatment, reduction = "umap", group.by = "stim")
p2 <- DimPlot(treatment, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(treatment, reduction = "umap", split.by = "stim", label = TRUE)
p1 + p2


#### Identify Conserved cell type markers #####
# performing differential expression after integration :
DefaultAssay(treatment) <- "RNA"
Cluster0Markers <- FindConservedMarkers(treatment, ident.1= 0, grouping.var = "stim", verbose = FALSE)
Cluster1Markers <- FindConservedMarkers(treatment, ident.1= 1, grouping.var = "stim", verbose = FALSE)
Cluster2Markers <- FindConservedMarkers(treatment, ident.1= 2, grouping.var = "stim", verbose = FALSE)
Cluster3Markers <- FindConservedMarkers(treatment, ident.1= 3, grouping.var = "stim", verbose = FALSE)



# output = list of genes 

# check for validity of markers within this list of genes : 
# using the markers identified, a set of Feature plots will br made fo each gene (marker) to see if they acurately describe clusters: 
Cluster0FP <-  FeaturePlot(treatment, features = c("CCR7", "CD7"), min.cutoff = "q9")
Cluster1FP <-  FeaturePlot(treatment, features = c("CCL2"))
Cluster2FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster3FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster4FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster5FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster6FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster7FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster8FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster9FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster10FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")
Cluster11FP <-  FeaturePlot(treatment, features = c(""), min.cutoff = "q9")

VlnPlot(treatment, features = "FCGR3A")

# After research, rename the clusters according to immune cell type : 

treatment <- RenameIdents(treatment, 
                          '0' = '',
                          '1' = '',
                          '2' = '',
                          '3' = '',
                          '4' = '',
                          '5' = '',
                          '6' = '',
                          '7' = '',
                          '8' = '',
                          '9' = '',
                          '10' = '',
                          '11' = '')

p4 <- DimPlot(treatment, label = TRUE)



# view conserved cell markers across conditions : 
# expression level + pergentage of cells in a cluster 

Idents(treatment) <- factor(Idents(treatment), 
                            levels = c("identified cell types", "")
                            )
Markers <- c("gene markers to plot", "") # top 3 marker genes for each cluster 

p5 <- DotPlot(treatment, 
              features = c("CD14", "CD68", "FCGR3A", "CCR2"), 
              cols = c("#41B6C4","#FF9AA2" ), 
              dot.scale = 8, 
              split.by = "stim")  + 
  RotatedAxis()

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
write.csv(predictions$labels,"honours/results/SRPredictionsBEDFine")
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
