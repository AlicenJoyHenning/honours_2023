##### Loading dependencies #####

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)

getwd()

alpha.data <-  Read10X(data.dir = "honours/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha.data, project='ifnalpha', min.cells=3, min.features=200)
# output : Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

alpha
# An object of class Seurat 
# 17692 features across 6158 samples within 1 assay 
# Active assay: RNA (17692 features, 0 variable features)

alpha.data[1:3, 1:3]
# 3 x 3 sparse Matrix of class "dgCMatrix"
# AAACCCAAGACTCATC AAACCCACAACGTATC AAACCCACAATTGCTG
# HGNC                  .                .                .
# PRDM16                .                .                .
# PEX10                 .                .                .

##### Removing unwanted cells based on # genes expressed and mito chondrial gene xpression #####

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")

alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

##### Normalizing the data #####

# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|

# 2 : feature selection 
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
# Calculating gene variances
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Calculating feature variances of standardized and clipped values
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|

# 3 : scaling the data 
all.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.genes)
# Centering and scaling data matrix
# |==========================================================================| 100%
# The results of this are stored in pbmc[["RNA"]]@scale.data


##### Perform Linear Dimensional Reduction #####
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input

alpha <- RunPCA(alpha, features = VariableFeatures(object = alpha))
# PC_ 1 
# Positive:  RPL23A, RPL3, RPL13A, RPL26, RPS3A, RPL18A, RPS23, RPL7A, RPS18, RPL35 
# ENSG00000237550, RPSA, RPL10A, RPL5, RPL10, RPS4X, EEF1B2, RPL6, RPS7, RPS8 
# RPS5, RPL18, RPS20, RPL37A, RPS15A, RPL12, RPL24, RPL35A, RPS6, RPS3 
# Negative:  IL1RN, S100A8, FTH1, FFAR2, MNDA, ACSL1, APOBEC3A, TYROBP, NCF1C, NINJ1 
# GLUL, SNX10, IFITM3, FCER1G, SPI1, NCF1, TYMP, C15orf48, IER3, TLR2 
# CCL4, CLEC7A, CCL4L2, LST1, S100A11, NCF2, VNN2, SAMSN1, ISG15, BASP1 
# PC_ 2 
# Positive:  LEF1, LDHB, IL7R, BCL11B, IKZF1, TCF7, GPR155, ADTRP, ITGA6, CAMK4 
# GIMAP4, CMTM8, LEPROTL1, C12orf57, MAL, PRMT2, PDCD4, CD3G, CCR7, TNFRSF25 
# MYC, SCML1, LTB, STING1, GIMAP5, FHIT, FAM184A, GIMAP6, PASK, ABLIM1 
# Negative:  HLA-DPA1, HLA-DRA, CD74, MYOF, KYNU, HLA-DQB1, LILRB1, CD86, RGL1, HLA-DQA1 
# LGALS1, ADAM28, EPB41L3, CPVL, CYP1B1, OLR1, ANXA2, LACC1, EMP1, LY86 
# BANK1, MS4A1, TCF4, MS4A7, RIN2, CCL2, CCL7, SWAP70, ENG, FGD2 
# PC_ 3 
# Positive:  MS4A1, IGHM, BANK1, BCL11A, IGHD, IGKC, HLA-DQA1, HLA-DPA1, ADAM28, HLA-DRA 
# TCF4, IGKJ5, CCSER1, CD19, PLEKHG7, RAB30, BCL2, MEF2C, NAPSB, EAF2 
# HLA-DQB1, CCR7, IGLC2, CD79B, HLA-DPB1, STAP1, CD24, ACSM3, CD74, C7orf50 
# Negative:  NKG7, PRF1, GNLY, ANXA1, GZMB, CLIC3, SAMD3, MYBL1, METRNL, SYTL2 
# CD8A, FYN, SH2D2A, C1orf21, KLRD1, IL18RAP, GZMM, CD96, PYHIN1, PLAAT4 
# SYTL3, AUTS2, EOMES, KLRF1, IQGAP2, PIK3R1, LAG3, MAF, TNFRSF18, TRGC2 
# PC_ 4 
# Positive:  NKG7, PRF1, GNLY, GZMB, CLIC3, SAMD3, KLRD1, MS4A1, KLRF1, BANK1 
# IGHM, BCL11A, IL18RAP, MYBL1, CD38, C1orf21, IGHD, SH2D2A, EOMES, ADAM28 
# CD8A, IGKC, STAP1, TNFRSF18, SYTL2, CCSER1, PYHIN1, FGFBP2, GZMM, AUTS2 
# Negative:  MYOF, EPB41L3, RGL1, CYP1B1, CPVL, EMP1, OLR1, CCL2, CCL7, LACC1 
# CD86, CXCL11, ENG, MS4A7, PLA2G7, CTSB, P2RY6, LYZ, MSR1, LILRB4 
# TGFBI, LILRB1, LEF1, SRC, LGALS1, CST3, KYNU, RIN2, THBS1, RGS10 
# PC_ 5 
# Positive:  ENSG00000284874, NRGN, ACRBP, TSC22D1, DAB2, PTCRA, TUBA4A, CMTM5, MMD, TMEM40 
# C2orf88, ENKUR, CLU, PGRMC1, CTTN, H2BC11, MPIG6B, MTURN, SPARC, PDE5A 
# GTF3C2, MAP3K7CL, H2AJ, TNFSF4, CLDN5, BEX3, H2AC6, CTSA, SSX2IP, TSPAN33 
# Negative:  ISG15, MYOF, RGL1, LGALS1, S100A11, APOBEC3A, LILRB1, TYROBP, OLR1, MS4A7 
# LACC1, RIN2, LILRB4, CPVL, KYNU, EMP1, IL1RN, ANXA5, CTSB, S100A4 
# IFITM3, ENG, CYP1B1, MT2A, CCL4, PLA2G7, S100A6, GCH1, IL4I1, CD86 


# Examine and visualize PCA results a few different ways :

print(alpha[["pca"]], dims = 1:5, nfeatures = 5)
# PC_ 1 
# Positive:  RPL23A, RPL3, RPL13A, RPL26, RPS3A 
# Negative:  IL1RN, S100A8, FTH1, FFAR2, MNDA 
# PC_ 2 
# Positive:  LEF1, LDHB, IL7R, BCL11B, IKZF1 
# Negative:  HLA-DPA1, HLA-DRA, CD74, MYOF, KYNU 
# PC_ 3 
# Positive:  MS4A1, IGHM, BANK1, BCL11A, IGHD 
# Negative:  NKG7, PRF1, GNLY, ANXA1, GZMB 
# PC_ 4 
# Positive:  NKG7, PRF1, GNLY, GZMB, CLIC3 
# Negative:  MYOF, EPB41L3, RGL1, CYP1B1, CPVL 
# PC_ 5 
# Positive:  ENSG00000284874, NRGN, ACRBP, TSC22D1, DAB2 
# Negative:  ISG15, MYOF, RGL1, LGALS1, S100A11 

DimHeatmap(alpha, dims = 1, cells = 500, balanced = TRUE)


DimPlot(alpha, reduction = "pca")

##### Determine dimensionality of the dataset ####

alpha <- JackStraw(alpha, num.replicate = 100)

alpha <- ScoreJackStraw(alpha, dims = 1:20)

JackStrawPlot(alpha, dims = 1:50)

##### Cluster cells ##### 

alpha <- FindNeighbors(alpha, dims = 1:20)
# Computing nearest neighbor graph
# Computing SNN
# the first 20 principal components are being used to calculate the nearest neighbors for the cells in the alpha dataset.

alpha.r25 <- FindClusters(alpha, resolution = 0.25) 
alpha <- FindClusters(alpha, resolution = 0.5) 
alpha.r65 <- FindClusters(alpha, resolution = 0.65) 
# clustering cells based on their similarity in a lower-dimensional space, typically after performing dimensionality reduction and finding nearest neighbors
# the cluster information (assignments) obtained from the FindClusters function is stored in the seurat object as a new metadata column. Each cell in your dataset is assigned to a specific cluster, and this assignment is added as a metadata attribute to the seurat object.

alpha.clusters <- as.data.frame(alpha@meta.data)



# /////

# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 4742
# Number of edges: 187903
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.8603
#    Number of communities: 11
#    Elapsed time: 0 seconds
# Look at cluster IDs of the first 5 cells
head(Idents(alpha), 5)
# AAACCCAAGACTCATC AAACCCACAATTGCTG AAACCCACAGCAAGAC AAACCCATCATTCATC 
# 1                3                1                2 
# AAACGAACAAAGGCTG 
# 2 
# Levels: 0 1 2 3 4 5 6 7 8 9 10


##### Run non-linear dimensional reduction #####
# using resolution 0.5 as my standard : 
alpha <- RunUMAP(alpha, dims = 1:10)
DimPlot(alpha, reduction = "umap")

alpha.r25 <- RunUMAP(alpha.r25, dims = 1:20)
DimPlot(alpha.r25, reduction = "umap")
# saveRDS(alpha.r25, file = "git_backup/plots/alpha.r25_UMAP")- too large to save 


alpha.r65 <- RunUMAP(alpha.r65, dims = 1:10)
DimPlot(alpha.r65, reduction = "umap")


# for the DimPlot, the clusters that are close together are presented in the same colour. This makes it difficult to view. To change this, 
# I want to se the R Color Brewer palette : 
install.packages("RColorBrewer")
library(RColorBrewer)
palette.a <- brewer.pal(11, "Paired")

Seurat::DimPlot(
  object = alpha,
  reduction = 'umap',
  group.by = 'seurat_clusters',
  pt.size = 1,
  label = FALSE,
  cols = palette.a
)


##### Finding differentially expressed features #####

cluster1.markers <- FindMarkers(alpha, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# S100A8   1.089219e-275   1.949636 0.975 0.437 1.927046e-271
# SLC25A37 5.253270e-223   1.562825 0.960 0.456 9.294084e-219
# LRRK2    4.276637e-217   1.673495 0.894 0.361 7.566227e-213
# NAMPTP1  3.880641e-207   1.444335 0.982 0.518 6.865630e-203
# VNN2     3.834474e-206   2.122229 0.610 0.152 6.783952e-202

# distinguishing between cluster 0 and 1 : 
cluster0v1.markers <- FindMarkers(alpha, ident.1 = 0, ident.2 = 1, min.pct =0.25)
head(cluster0v1.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# VNN2   4.176859e-59 -1.2197463 0.295 0.610 7.389698e-55
# S100A8 2.409076e-56 -0.7544300 0.898 0.975 4.262137e-52
# CCL4L2 1.115124e-42  1.3563812 0.566 0.297 1.972877e-38
# CCL4   1.665484e-40  1.1984775 0.665 0.418 2.946575e-36
# ACTB   4.098977e-40 -0.7574825 0.675 0.862 7.251911e-36

cluster0.markers <- FindMarkers(alpha, ident.1 = 0, ident.2 = c(1, 10), min.pct = 0.25)
head(cluster0.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# VNN2   1.286096e-54 -1.1800757 0.295 0.593 2.275361e-50
# S100A8 8.568999e-47 -0.7103968 0.898 0.951 1.516027e-42
# CCL4L2 7.527496e-45  1.3885013 0.566 0.295 1.331765e-40
# ACTB   1.014142e-44 -0.7992654 0.675 0.866 1.794220e-40
# CCL4   3.556346e-43  1.2325919 0.665 0.413 6.291888e-39

# find markers for every cluster compared to all remaining cells, report only the positive nes
alpha.r65.markers <- FindAllMarkers(alpha.r65, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

alpha.r65.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

# 0.5 (round 1) : CCL4L2, CCL4, VNN2, S100A8,  CCR7, RPS13 , MAF,  IL32 ,  PASK, NPM1 
# 0.25 : IL1RN, SOD2, CCR7, RPS13, MAF, IL32, PASK, NPM1, IGKC, CD74, GNLY, GZMB, PRF1, CD8A, CCL2, CCL7, ATP10D, HDC, NRGN
# 0.5 : CCL4L2,  CCL4,  VNN2, S100A8, CCR7, RPS13,  MAF, IL32, PASK,NPM1  , 
# 0.65 : CCL4L2, CCL4, IDO1, GBP5 , CCR7, RPS13, MAF, IL32,VNN2, S100A8   

# using the markers identified, a set of Feature plots will br made fo each gene (marker) to see if they acurately describe clusters: 
FeaturePlot(alpha.r5, features = c("LYZ","FCER1A"))

##### Create data frame that stores all information about the differentially expressed features of each cluster ##### 

cluster = c(seq(0, 10)) # column 1
unique.markers = FindMarkers(alpha, ident.1 = i, min.pct = 0.25) | head n = 5
cluster.comparision.markers <- FindMarkers(alpha, ident.1 = i, ident.2 = j, min.pct =0.25)

diff.features <- data.frame(cluster, 
                            unique.markers, 
                            cluster.comparision.markers)


##### Annotate clusters #####

# Using the scAnnotatR package with pretrained models to classify cell types 

BiocManager::install("devtools")
BiocManager::install("scAnnotatR")

library(scAnnotatR)

default_models <- load_models("default") # loading the pre-trained models 
names(default_models)
# [1] "B cells"           "Plasma cells"      "NK"                "CD16 NK"           "CD56 NK"           "T cells"          
# [7] "CD4 T cells"       "CD8 T cells"       "Treg"              "NKT"               "ILC"               "Monocytes"        
# [13] "CD14 Mono"         "CD16 Mono"         "DC"                "pDC"               "Endothelial cells" "LEC"              
# [19] "VEC"               "Platelets"         "RBC"               "Melanocyte"        "Schwann cells"     "Pericytes"        
# [25] "Mast cells"        "Keratinocytes"     "alpha"             "beta"              "delta"             "gamma"            
# [31] "acinar"            "ductal"            "Fibroblasts"

# takes as input a seurat object : 
is(alpha.r5, "Seurat")
# [1] TRUE


# To launch cell type identification, we simply call the `classify_cells`function : 
alpha.r5.scannotatR <- classify_cells(classify_obj = alpha.r5, 
                             assay = 'RNA', slot = 'counts',
                             cell_types = c('NK', ), 
                             path_to_models = 'default')

#  returns the input object but with additional columns in the metadata table.
# New columns are:
#   
#   * **predicted_cell_type**: The predicted cell type, also containing any 
# ambiguous assignments. In these cases, the possible cell types are separated
# by a "/"
# 
# * **most_probable_cell_type**: contains the most probably cell type ignoring any 
# ambiguous assignments.
# 
# * columns with syntax `[celltype]_p`: probability of a cell to belong 
# to a cell type. Unknown cell types are marked as NAs.

# The predicted cell types can now simply be visualized using the matching plotting functions

Seurat::DimPlot(alpha.r5.scannotatR, group.by = "most_probable_cell_type")

# For a certain cell type, users can also view the prediction probability

Seurat::FeaturePlot(seurat.obj, features = "B_cells_p")








