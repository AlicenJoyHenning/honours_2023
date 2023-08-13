# Dataset Processing in R with Seurat 
https://satijalab.org/seurat/articles/pbmc3k_tutorial


## Preprocessing Stages  
### _Loading the datasets_

After going through the first two functions of SASCRiP, the output folder containing barcodes, features, and matrix files were used to go through a manual preprocessing workflow offered by Seurat. 

The ```Read10X()``` function from 10X Genomics takes as input the ouptput of 10X cellranger and returns a UMI count matrix: a matrix that represents the number of molecules for each feature (meaning gene > row) that is detected in each cell (column). The count matrix is then used to create a ```Seurat object``` that stores raw and analysis data for a single cell dataset.

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

### _Quality control_
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
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/065e352d-7f02-4171-90dd-fd81b050ce93)

### _Normalization_

After quality control, the dataset must be normalized. 

**1) Log Normalization**

This means applying a logarithm transformation to the raw counts of gene expression to reduce the impact of high-variance genes. High-variance genes are those that exhibit large differences in expression levels across cells. In scRNA-seq datasets, some genes might have very high expression in a few cells while being expressed at lower levels in others. These high-variance genes can introduce noise and dominate the overall variability in the dataset. Log-normalization helps reduce the impact of high-variance genes by compressing their expression values. Since the logarithm function compresses high values more than low values, the resulting log-transformed values tend to have a more balanced distribution, which reduces the dominance of extreme values. Log normalization also reduces the effect of differences in library sizes (total counts) across cells. 

Log-normalization will be applied using a function in the Seurat package (LogNormalize). Log-transforming the data allows for the representation of fold changes in expression on a linear scale. A 2-fold change in expression corresponds to a difference of 1 unit on the logarithmic scale. This is particularly useful for visualizations and downstream analyses. A 10 000 scaling factor (default) was chosen. 

**2) Finding variable genes**

The functions FindVariableFeatures() with the selection.method = "vst" argument will identify e a subset of highly variable genes (features) from the RNA matrix data stored in the Seurat object. These functions will result in a reduced number of genes/features for downstream analyses - it it does not alter the total number of genes, ```dim(object@assays$RNA@counts0[2]```, in your original expression matrix, only adds new information to the variable gene identification : ```length(object@assays$RNA@var.features)```

The function FindVariableFeatures() from the Seurat package in R is used to identify highly variable features (genes) in a dataset, essentially identifying the most informative genes that contribute to cell-to-cell variability. . Highly variable features is important because not all genes have the same level of biological variability across cells. Some genes are highly expressed and exhibit significant variation, while others have lower expression and variability. Identifying these highly variable features helps focus downstream analyses on genes that carry more biologically relevant information.

This used the "vst" (Variance Stabilizing Transformation), which calculates the coefficient of variation (CV) and uses it to measure variability in gene expression. The selection.method = "vst" argument indicates that the function will calculate the coefficient of variation (CV) after applying a variance stabilizing transformation (VST) to the data. The VST helps stabilize the variance across the dynamic range of expression, making the data more suitable for detecting true biological variability. nfeatures sets the maximum number of variable features to identify. In this case, it's set to 2000, indicating that the function will identify up to 2000 most variable features.

 **3) Scaling the data**

The function ScaleData() from the Seurat package in R is used to perform data scaling on the gene expression data. Scaling the data involves centering and rescaling each gene's expression values in a way that makes them comparable across cells. This is important because the absolute expression values can vary widely between cells due to technical factors such as sequencing depth or capture efficiency. Scaling ensures that the differences in expression are driven by biology rather than technical variability. This transformation centers the data around the mean of each gene, and it rescales the data to have a consistent unit of measurement (standard deviation).

Scaling data is often a prerequisite for downstream analyses like principal component analysis (PCA), t-distributed stochastic neighbor embedding (t-SNE), and clustering. These analyses benefit from scaled data because they focus on patterns of relative expression rather than absolute expression levels.

![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/80e02fda-0b34-4c62-80ae-a126e650861d)


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

### _Linear dimensional reduction_
## Linear dimensional reduction 

Next perform PCA on the scaled data. By default, only the previously determined variable features (```object@assays$RNA@var.features```) are used as input. The function RunPCA() from the Seurat package in R is used to perform principal component analysis (PCA) on the gene expression data.  PCA reduces the dimensionality of the gene expression data while retaining the most informative patterns of variation. It transforms the original high-dimensional gene expression space into a new set of orthogonal (uncorrelated) axes called principal components. Each principal component captures different sources of variability in the data.

By retaining only the top principal components that explain the most variance, PCA helps to identify key patterns such as biological variability, batch effects, or other sources of variation across cells. PCA is commonly used as a preprocessing step for downstream analyses like clustering, visualization (e.g., t-SNE), and differential gene expression analysis. It provides a reduced-dimensional representation of the data that is more manageable and informative for these analyses.

The features argument specifies which genes will be used to perform PCA. By setting features to the variable features identified using VariableFeatures(), you are focusing PCA on the most informative genes that exhibit high biological variability.

Performing PCA using the RunPCA() function in Seurat will add new attributes to the Seurat object and store the results of the PCA analysis. 


```R
alpha <- RunPCA(alpha, features = VariableFeatures(object = alpha))
lambda <- RunPCA(lambda, features = VariableFeatures(object = lambda))
untreated <- RunPCA(untreated, features = VariableFeatures(object = untreated))
```

This generates an output that gives the PC components (1:5) and lists underneath it the Positive and Negative genes (30 each). This is stored in the pca slot that can be viewed as follows, ```print(alpha[["pca"]], dims = 1:5, nfeatures = 5)```. Alternatively, you can manually enter the output as a dataframe stored in the seurat object under the assays section (see script : PCA_additions.R ) 

You can visualize the effect of PCA using various techniques such as scatter plots, t-SNE, or UMAP, which will help you see how cells are distributed in the reduced-dimensional space. Seurat provides several useful ways of visualizing both cells and features that define the PCA

_PCA Scatter Plots_
Each point represents a cell that is projected against the top principal components. 
```DimPlot(object, reduction ="pca", cols = c("blue"))```
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/3d458aac-ddc4-4eef-bf12-1206689978dd)

## Determine dimensionality of the dataset 

In scRNA-seq data, there can be a lot of noise (variability not related to biological differences) in the gene expression measurements. To make sense of the data and find meaningful patterns, it's important to reduce this noise. One way to do this is by using Principal Component Analysis (PCA), which finds the most important combinations of genes that capture the major sources of variation in the data.
After performing dimensionality reduction (like PCA) on your single-cell data, you end up with principal components that capture different levels of information. Some PCs might capture meaningful biological variation, while others might capture random noise. By running JackStraw, you can identify which PCs have more significant variation than expected by chance alone. 

The statistical test known as "JackStraw" can be performed on a seurat object that helps you determine if any of the principal components are significant when compared to random noise. After the statistical test has been performed, the _score_ for each principal component based on how much it deviates from random noise




```R
alpha <- JackStraw(alpha, num.replicate = 100)
alpha <- ScoreJackStraw(alpha, dims = 1:20)
js_p1 <- JackStrawPlot(alpha, dims = 1:20)


lambda <- JackStraw(lambda, num.replicate = 100)
lambda <- ScoreJackStraw(lambda, dims = 1:20)
js_p2 <- JackStrawPlot(lambda, dims = 1:20)

untreated <- JackStraw(untreated, num.replicate = 100)
untreated <- ScoreJackStraw(untreated, dims = 1:20)
js_p3 <- JackStrawPlot(untreated, dims = 1:20)

grid.arrange(js_p1, js_p2, js_p3, ncol = 3)
```
The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). JackStraw Procedure: The JackStraw procedure involves randomly shuffling the data a bit (permuting) and then re-running PCA. By doing this multiple times, you create a "null distribution" of how the scores might look if there's no real signal. PCs that significantly deviate from this null distribution are considered significant.

It shows the distribution of p-values for each PC. A p-value is a measure of how likely it is to observe a certain value by random chance. If a PC's p-value is very low (solid curve above the dashed line), it means that the observed score is unlikely to have occurred by random chance alone. These significant PCs are the ones you should pay attention to. Lower p-values indicate more significant PCs, suggesting that the genes driving those PCs have biological relevance. This plot helps you identify a "gap" where the observed p-values separate from random noise p-values, indicating the PCs that are likely to be biologically relevant.

By identifying the significant PCs, you're ensuring that the features (genes) you're using for clustering cells are meaningful and not dominated by noise. This step helps you focus on the most important sources of variation and makes your subsequent analyses, like clustering, more accurate and biologically relevant.



![image](


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




_t-SNE Plots_
Dimensionality reduction technique that highlights cells in close proximity (well not as good at capturing global structure as UMAP). 
```R
alpha <- RunTSNE(alpha)
lambda <- RunTSNE(lambda)
untreated <- RunTSNE(untreated)

tsne_p1 <- DimPlot(alpha, reduction = "tsne", cols = "darkblue")
tsne_p2 <- DimPlot(lambda, reduction = "tsne", cols = "blue")
tsne_p3 <- DimPlot(untreated, reduction = "tsne", cols = "lightblue")
library(gridExtra)
grid.arrange(tsne_p1, tsne_p2, tsne_p3, ncol = 3)
```
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/5b51811c-397e-400d-a35f-06e8cf8b00f5)

_UMAP Plots_
UMAP (Uniform Manifold Approximation and Projection) is another dimensionality reduction technique that often provides better separation and global structure preservation than t-SNE. UMAP is especially helpful for revealing subtle population differences and identifying clusters.
```R
alpha <- RunUMAP(alpha, dims = 1:10)
lambda <- RunUMAP(lambda, dims = 1:10)
untreated <- RunUMAP(untreated, dims = 1:10)

umap_p1 <- DimPlot(alpha, reduction = "umap", cols = "darkblue")
umap_p2 <- DimPlot(lambda, reduction = "umap", cols = "blue")
umap_p3 <- DimPlot(untreated, reduction = "umap", cols = "lightblue")
library(gridExtra)
grid.arrange(umap_p1, umap_p2, umap_p3, ncol = 3)
```
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/35744a8a-964e-4b43-870a-5691077975c3)

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
