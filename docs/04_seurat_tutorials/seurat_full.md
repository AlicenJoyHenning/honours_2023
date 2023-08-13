# Dataset Processing in R with Seurat 
https://satijalab.org/seurat/articles/pbmc3k_tutorial


## Preprocessing Stages  
After going through the first two functions of SASCRiP, the output folder containing barcodes, features, and matrix files were used to go through a manual preprocessing workflow offered by Seurat.

```R
# Load the dependencies :

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)

# Load the data :

alpha.data <-  Read10X(data.dir = "honours/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha.data, project='ifnalpha', min.cells=3, min.features=200)

lambda.data <-  Read10X(data.dir = "honours/ifnlambda/seurat_matrix/")
lambda <- CreateSeuratObject(counts=lambda.data, project='ifnlambda', min.cells=3, min.features=200)

untreated.data <-  Read10X(data.dir = "honours/untreated/seurat_matrix/")
untreated <- CreateSeuratObject(counts=untreated.data, project='untrearted', min.cells=3, min.features=200)

# output : Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
```

To view and explore the Seurat objects : 

```R 

alpha
# An object of class Seurat 
# 17692 features across 6158 samples within 1 assay 
# Active assay: RNA (17692 features, 0 variable features)

alpha.data[1:3, 1:3]
# 3 x 3 sparse Matrix of class "dgCMatrix"
# AAACCCAAGACTCATC AAACCCACAACGTATC AAACCCACAATTGCTG
# HGNC                  .                .                .
# PRDM16                .                .                .
# PEX10       
```

Next, the datasets must be altered to remove low quality cells, determined by the number of features expressed and the mitochondrial percentage. 

1. nFeature_RNA > 200:
This condition filters out cells that have very low RNA content. Cells with a low number of detected features (genes) might represent low-quality cells or technical artifacts.

2. nFeature_RNA < 2500:
This condition filters out cells with very high RNA content. Extremely high feature counts might indicate doublets (two cells captured together) or other technical issues.

3. percent.mt < 10:
This condition filters out cells with a high percentage of mitochondrial gene expression. High levels of mitochondrial gene expression can indicate poor cell health or cell stress. Cells with high mitochondrial content are often considered to be damaged or apoptotic.

```R
First, the mitochondrial genes are identified  by the Serat function : 

seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")
# The [[ operator adds columns to object metadata to stash QC stats
# PercentageFeatureSet() takes the alpha object as input, the pattern = "^MT-" argument specifies a regular expression pattern to identify mitochondrial genes. The "^MT-" pattern will match gene names that start with "MT-" (naming convention for mitochondrial genes)

# Adding metadata column :

alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-") # 6483 cells
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-") # 5922

# Doing the subsetting (quality control) : (checking the change in dim(object@assays$RNA@counts)[1]/[2]

alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # 17692 4742
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # 
18272 4811
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10) # 18188 5349
```
![image]("git_backup/images/cell_quality_control.jpg")


After quality control, the dataset must be normalized. 

** 1) Log Normalization**

This means applying a logarithm transformation to the raw counts of gene expression to reduce the impact of high-variance genes. High-variance genes are those that exhibit large differences in expression levels across cells. In scRNA-seq datasets, some genes might have very high expression in a few cells while being expressed at lower levels in others. These high-variance genes can introduce noise and dominate the overall variability in the dataset. Log-normalization helps reduce the impact of high-variance genes by compressing their expression values. Since the logarithm function compresses high values more than low values, the resulting log-transformed values tend to have a more balanced distribution, which reduces the dominance of extreme values. Log normalization also reduces the effect of differences in library sizes (total counts) across cells. 

Log-normalization will be applied using a function in the Seurat package (LogNormalize). Log-transforming the data allows for the representation of fold changes in expression on a linear scale. A 2-fold change in expression corresponds to a difference of 1 unit on the logarithmic scale. This is particularly useful for visualizations and downstream analyses. A 10 000 scaling factor (default) was chosen. 

** 2) Finding variable genes**

The functions FindVariableFeatures() with the selection.method = "vst" argument will identify e a subset of highly variable genes (features) from the RNA matrix data stored in the Seurat object. These functions will result in a reduced number of genes/features for downstream analyses - it it does not alter the total number of genes, ```dim(object@assays$RNA@counts0[2]```, in your original expression matrix, only adds new information to the variable gene identification : ```length(object@assays$RNA@var.features)```

The function FindVariableFeatures() from the Seurat package in R is used to identify highly variable features (genes) in a dataset, essentially identifying the most informative genes that contribute to cell-to-cell variability. . Highly variable features is important because not all genes have the same level of biological variability across cells. Some genes are highly expressed and exhibit significant variation, while others have lower expression and variability. Identifying these highly variable features helps focus downstream analyses on genes that carry more biologically relevant information.

This used the "vst" (Variance Stabilizing Transformation), which calculates the coefficient of variation (CV) and uses it to measure variability in gene expression. The selection.method = "vst" argument indicates that the function will calculate the coefficient of variation (CV) after applying a variance stabilizing transformation (VST) to the data. The VST helps stabilize the variance across the dynamic range of expression, making the data more suitable for detecting true biological variability. nfeatures sets the maximum number of variable features to identify. In this case, it's set to 2000, indicating that the function will identify up to 2000 most variable features.

 ** 3) Scaling the data **

The function ScaleData() from the Seurat package in R is used to perform data scaling on the gene expression data. Scaling the data involves centering and rescaling each gene's expression values in a way that makes them comparable across cells. This is important because the absolute expression values can vary widely between cells due to technical factors such as sequencing depth or capture efficiency. Scaling ensures that the differences in expression are driven by biology rather than technical variability. This transformation centers the data around the mean of each gene, and it rescales the data to have a consistent unit of measurement (standard deviation).

Scaling data is often a prerequisite for downstream analyses like principal component analysis (PCA), t-distributed stochastic neighbor embedding (t-SNE), and clustering. These analyses benefit from scaled data because they focus on patterns of relative expression rather than absolute expression levels.

![image]("")

```R

# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
lambda <- NormalizeData(lambda, normalization.method = "LogNormalize", scale.factor = 10000)
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)


# 2 : feature selection 
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000) #
lambda <- FindVariableFeatures(lambda, selection.method = "vst", nfeatures = 2000)
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)

# 3 : scaling the data 
all.alpha.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.alpha.genes)
all.lambda.genes <- rownames(lambda)
lambda <- ScaleData(lambda, features = all.lambda.genes)
all.untreated.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.untreated.genes)

```



## Linear dimensional reduction 

Next perform PCA on the scaled data. By default, only the previously determined variable features are used as input

```R alpha <- RunPCA(alpha, features = VariableFeatures(object = alpha))```

This generates the output: 
```R
PC_ 1 
Positive:  RPL23A, RPL3, RPL13A, RPL26, RPS3A, RPL18A, RPS23, RPL7A, RPS18, RPL35 
ENSG00000237550, RPSA, RPL10A, RPL5, RPL10, RPS4X, EEF1B2, RPL6, RPS7, RPS8 
RPS5, RPL18, RPS20, RPL37A, RPS15A, RPL12, RPL24, RPL35A, RPS6, RPS3 
Negative:  IL1RN, S100A8, FTH1, FFAR2, MNDA, ACSL1, APOBEC3A, TYROBP, NCF1C, NINJ1 
GLUL, SNX10, IFITM3, FCER1G, SPI1, NCF1, TYMP, C15orf48, IER3, TLR2 
CCL4, CLEC7A, CCL4L2, LST1, S100A11, NCF2, VNN2, SAMSN1, ISG15, BASP1 
PC_ 2 
Positive:  LEF1, LDHB, IL7R, BCL11B, IKZF1, TCF7, GPR155, ADTRP, ITGA6, CAMK4 
GIMAP4, CMTM8, LEPROTL1, C12orf57, MAL, PRMT2, PDCD4, CD3G, CCR7, TNFRSF25 
MYC, SCML1, LTB, STING1, GIMAP5, FHIT, FAM184A, GIMAP6, PASK, ABLIM1 
Negative:  HLA-DPA1, HLA-DRA, CD74, MYOF, KYNU, HLA-DQB1, LILRB1, CD86, RGL1, HLA-DQA1 
LGALS1, ADAM28, EPB41L3, CPVL, CYP1B1, OLR1, ANXA2, LACC1, EMP1, LY86 
BANK1, MS4A1, TCF4, MS4A7, RIN2, CCL2, CCL7, SWAP70, ENG, FGD2 
PC_ 3 
Positive:  MS4A1, IGHM, BANK1, BCL11A, IGHD, IGKC, HLA-DQA1, HLA-DPA1, ADAM28, HLA-DRA 
TCF4, IGKJ5, CCSER1, CD19, PLEKHG7, RAB30, BCL2, MEF2C, NAPSB, EAF2 
HLA-DQB1, CCR7, IGLC2, CD79B, HLA-DPB1, STAP1, CD24, ACSM3, CD74, C7orf50 
Negative:  NKG7, PRF1, GNLY, ANXA1, GZMB, CLIC3, SAMD3, MYBL1, METRNL, SYTL2 
CD8A, FYN, SH2D2A, C1orf21, KLRD1, IL18RAP, GZMM, CD96, PYHIN1, PLAAT4 
SYTL3, AUTS2, EOMES, KLRF1, IQGAP2, PIK3R1, LAG3, MAF, TNFRSF18, TRGC2 
PC_ 4 
Positive:  NKG7, PRF1, GNLY, GZMB, CLIC3, SAMD3, KLRD1, MS4A1, KLRF1, BANK1 
IGHM, BCL11A, IL18RAP, MYBL1, CD38, C1orf21, IGHD, SH2D2A, EOMES, ADAM28 
CD8A, IGKC, STAP1, TNFRSF18, SYTL2, CCSER1, PYHIN1, FGFBP2, GZMM, AUTS2 
Negative:  MYOF, EPB41L3, RGL1, CYP1B1, CPVL, EMP1, OLR1, CCL2, CCL7, LACC1 
CD86, CXCL11, ENG, MS4A7, PLA2G7, CTSB, P2RY6, LYZ, MSR1, LILRB4 
TGFBI, LILRB1, LEF1, SRC, LGALS1, CST3, KYNU, RIN2, THBS1, RGS10 
PC_ 5 
Positive:  ENSG00000284874, NRGN, ACRBP, TSC22D1, DAB2, PTCRA, TUBA4A, CMTM5, MMD, TMEM40 
C2orf88, ENKUR, CLU, PGRMC1, CTTN, H2BC11, MPIG6B, MTURN, SPARC, PDE5A 
GTF3C2, MAP3K7CL, H2AJ, TNFSF4, CLDN5, BEX3, H2AC6, CTSA, SSX2IP, TSPAN33 
Negative:  ISG15, MYOF, RGL1, LGALS1, S100A11, APOBEC3A, LILRB1, TYROBP, OLR1, MS4A7 
LACC1, RIN2, LILRB4, CPVL, KYNU, EMP1, IL1RN, ANXA5, CTSB, S100A4 
IFITM3, ENG, CYP1B1, MT2A, CCL4, PLA2G7, S100A6, GCH1, IL4I1, CD86 
```
This code performs Principal Component Analysis (PCA) on a Seurat object named alpha.
alpha: This is the Seurat object on which you want to perform PCA.
RunPCA(): This is a function provided by the Seurat package that performs PCA. PCA is a dimensionality reduction technique that is commonly used to reduce the complexity of high-dimensional data while retaining as much information as possible.
features = VariableFeatures(object = alpha): This part specifies the features (genes) to be used for PCA. It uses the VariableFeatures() function to select the highly variable features from the Seurat object alpha. Highly variable features are those that exhibit the most biological variation across cells, making them important for capturing cell-to-cell differences.

Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. 

![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/9f67f1ab-1766-48cf-b795-b68ee890a7b6)

## Determine dimensionality of the dataset 

In scRNA-seq data, there can be a lot of noise (variability not related to biological differences) in the gene expression measurements. To make sense of the data and find meaningful patterns, it's important to reduce this noise. One way to do this is by using Principal Component Analysis (PCA), which finds the most important combinations of genes that capture the major sources of variation in the data.

```R
alpha <- JackStraw(alpha, num.replicate = 100)
alpha <- ScoreJackStraw(alpha, dims = 1:20)
```
The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). JackStraw Procedure: The JackStraw procedure involves randomly shuffling the data a bit (permuting) and then re-running PCA. By doing this multiple times, you create a "null distribution" of how the scores might look if there's no real signal. PCs that significantly deviate from this null distribution are considered significant.

It shows the distribution of p-values for each PC. A p-value is a measure of how likely it is to observe a certain value by random chance. If a PC's p-value is very low (solid curve above the dashed line), it means that the observed score is unlikely to have occurred by random chance alone. These significant PCs are the ones you should pay attention to.

By identifying the significant PCs, you're ensuring that the features (genes) you're using for clustering cells are meaningful and not dominated by noise. This step helps you focus on the most important sources of variation and makes your subsequent analyses, like clustering, more accurate and biologically relevant.

```R
JackStrawPlot(alpha, dims = 1:20)
```
![image](https://github.com/AlicenJoyHenning/honours_2023/blob/main/plots/alpha_JackStraw.jpg)
Using maximum dimensions it seems like all the PCAs are significant?

## Cluster cells 

Seurat v3 applies a graph-based clustering approach, building upon initial strategies in (Macosko et al). Importantly, the distance metric which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partitioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [SNN-Cliq, Xu and Su, Bioinformatics, 2015] and CyTOF data [PhenoGraph, Levine et al., Cell, 2015]. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors() function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).
```R
alpha <- FindNeighbors(alpha, dims = 1:20)
# Computing nearest neighbor graph
# Computing SNN
```

To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters() function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. 

```R
alpha <- FindClusters(alpha, resolution = 0.5) 
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
```
Experimenting with different resolutions: 
In the context of single-cell RNA sequencing (scRNA-seq) data analysis and clustering using community detection algorithms, such as the Louvain method or Leiden algorithm, the concept of "resolution" has a specific meaning that can sometimes be a bit counterintuitive.

In these algorithms, the resolution parameter is used to control the granularity of the clusters that are identified in the data. A lower resolution value results in larger, more inclusive clusters, while a higher resolution value leads to smaller, more distinct clusters.


## Run non-linear dimensional reduction (UMAP/tSNE) 

Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

```R
alpha <- RunUMAP(alpha, dims = 1:20)
DimPlot(alpha, reduction = "umap")
saveRDS(alpha, file = "Git/plots/alpha_UMAP.rds")
```
![image](https://github.com/AlicenJoyHenning/honours_2023/blob/main/plots/alpha_UMAP_0.5.jpg?raw=true)

The colours and aesthetic of the DimPlot are not easily manipulated, so ggplot can be used as an alternative as follows : 

```R
# Alternatively, using ggplot :   

# Extract UMAP coordinates and cluster information
alpha.umap.coords <- as.data.frame(alpha@reductions$umap@cell.embeddings)
clusters <- alpha$seurat_clusters

# Create a dataframe for ggplot
alpha.df <- data.frame(
  x = alpha.umap.coords$UMAP_1,
  y = alpha.umap.coords$UMAP_2,
  seurat_clusters = factor(clusters)
)

# Define color palette
palette.a <- RColorBrewer::brewer.pal(11, "Paired")

# Create the ggplot plot
ggplot(alpha.df, aes(x, y, colour = seurat_clusters)) +
  geom_point(size = 1) +
  scale_colour_manual(values = palette.a) +
  labs(#title = "IFN alpha",
       x = "UMAP 1",  # Rename x-axis label
       y = "UMAP 2",
       color = "")  + 
  theme_minimal() + 
  theme(#panel.background = element_rect(fill = "lightgrey"),  # Set background color
       panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  # Increase axis label size
        axis.title = element_text(size = 14), 
        plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
        panel.border = element_rect(color = "black", fill = NA),
        #legend.background = element_rect(color = "black", fill = "white"),
        legend.position = "right", 
        legend.title = element_text(size =  14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))) + 
  guides(color = guide_legend(
    override.aes = list(
      #hape = rep(22, length(palette.a)),  # Use squares (blocks)
      fill = palette.a, 
      size = 3.5),  # Color the squares with the same palette
      key_height = unit(1, "npc"),  # Spread the legend dots across the vertical length
      key_width = unit(4, "cm"),   # Adjust the width of the legend blocks
      title.theme = element_text(hjust = 0.5),  # Center the legend title
    label.position = "right",
    label.hjust = 1
  ))
```
![image](https://github.com/AlicenJoyHenning/honours_2023/blob/main/plots/alpha.UMAP.good.colours.jpg)


## Finding differentially expressed features (cluster biomarkers) 
Seurat can help you find markers that define clusters via differential expression. By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

Using the  ```R FindMarkers()``` function, you can find the top markers for each cluster : 
```R
cluster1.markers <- FindMarkers(alpha, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2     p_val_adj
# S100A8   1.089219e-275   1.949636 0.975 0.437 1.927046e-271
# SLC25A37 5.253270e-223   1.562825 0.960 0.456 9.294084e-219
# LRRK2    4.276637e-217   1.673495 0.894 0.361 7.566227e-213
# NAMPTP1  3.880641e-207   1.444335 0.982 0.518 6.865630e-203
# VNN2     3.834474e-206   2.122229 0.610 0.152 6.783952e-202
```

Alternatively, you can look at the markers that distinguish one cluster from another. Here, we are distinguishing cluster 0 from cluster 1 (as they have grouped together on the plot) 

```R
# Distinguishing between cluster 0 and 1 : 
cluster0v1.markers <- FindMarkers(alpha, ident.1 = 0, ident.2 = 1, min.pct =0.25)
head(cluster0v1.markers, n = 5)

# p_val avg_log2FC pct.1 pct.2    p_val_adj
# VNN2   4.176859e-59 -1.2197463 0.295 0.610 7.389698e-55
# S100A8 2.409076e-56 -0.7544300 0.898 0.975 4.262137e-52
# CCL4L2 1.115124e-42  1.3563812 0.566 0.297 1.972877e-38
# CCL4   1.665484e-40  1.1984775 0.665 0.418 2.946575e-36
# ACTB   4.098977e-40 -0.7574825 0.675 0.862 7.251911e-36

# This can also be done to distinguish one cluster from a number of other clusters :

cluster0.markers <- FindMarkers(alpha, ident.1 = 0, ident.2 = c(1, 10), min.pct = 0.25)
head(cluster0.markers, n = 5)
# p_val avg_log2FC pct.1 pct.2    p_val_adj
# VNN2   1.286096e-54 -1.1800757 0.295 0.593 2.275361e-50
# S100A8 8.568999e-47 -0.7103968 0.898 0.951 1.516027e-42
# CCL4L2 7.527496e-45  1.3885013 0.566 0.295 1.331765e-40
# ACTB   1.014142e-44 -0.7992654 0.675 0.866 1.794220e-40
# CCL4   3.556346e-43  1.2325919 0.665 0.413 6.291888e-39


# find markers for every cluster compared to all remaining cells, report only the positive nes
alpha.markers <- FindAllMarkers(alpha, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

alpha.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
# CCL4L2, CCL4, VNN2, S100A8,  CCR7, RPS13 , MAF,  IL32 ,  PASK, NPM1


# using the markers identified, a set of Feature plots will br made fo each gene (marker) to see if they acurately describe clusters: 
FeaturePlot(alpha, features = c("CCL4L2", "CCL4", "VNN2", "S100A8", "CCR7", "RPS13", "MAF", "IL32", "PASK", "NPM1"))

# Alternatively, we can use DoHeatMap to show the top 20 markers for each cluster: (not working for me)

 top10 <- alpha.markers   %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

![image](https://github.com/AlicenJoyHenning/honours_2023/blob/main/plots/alpha_featuresplot_1.jpg)
