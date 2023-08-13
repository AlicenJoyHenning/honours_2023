##### Loading dependencies #####

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)

alpha.data <-  Read10X(data.dir = "honours/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha.data, project='ifnalpha', min.cells=3, min.features=200)
# output : Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

lambda.data <-  Read10X(data.dir = "honours/ifnlambda/seurat_matrix/")
lambda <- CreateSeuratObject(counts=lambda.data, project='ifnlambda', min.cells=3, min.features=200)

untreated.data <-  Read10X(data.dir = "honours/untreated/seurat_matrix/")
untreated <- CreateSeuratObject(counts=untreated.data, project='untrearted', min.cells=3, min.features=200)

# To save the seurat objects :  

saveRDS(alpha, file= "honours/ifnalpha/alpha.rds")
saveRDS(lambda, file= "honours/ifnlambda/lambda.rds")
saveRDS(untreated, file= "honours/untreated/untreated.rds")

# To load the saved seurat objects : 

alpha <- readRDS("honours/ifnalpha/alpha.rds")
lambda <- readRDS("honours/ifnlambda/lambda.rds")
untreated <- readRDS("honours/untreated/untreated.rds")

##### Removing unwanted cells based on # genes expressed and mito chondrial gene xpression #####


alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
lambda <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)


##### Normalizing the data #####

# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
lambda <- NormalizeData(lambda, normalization.method = "LogNormalize", scale.factor = 10000)
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)


# 2 : feature selection (check dim(object@assays$RNA@counts)[2] to view genes and length(oject@assays$RNA@var.features) to view 2000 variable genes)
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
lambda <- FindVariableFeatures(lambda, selection.method = "vst", nfeatures = 2000)
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)

# 3 : scaling the data 
# The results of this are stored in object[["RNA"]]@scale.data

all.alpha.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.alpha.genes)
all.lambda.genes <- rownames(lambda)
lambda <- ScaleData(lambda, features = all.lambda.genes)
all.untreated.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.untreated.genes)


##### Perform Linear Dimensional Reduction #####
# Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input

alpha <- RunPCA(alpha, features = VariableFeatures(object = alpha))
lambda <- RunPCA(lambda, features = VariableFeatures(object = lambda))
untreated <- RunPCA(untreated, features = VariableFeatures(object = untreated))

# Examine and visualize PCA results a few different ways 
DimHeatmap(alpha, dims = 1, cells = 500, balanced = TRUE)

# scatter plot : 
pca_p1 <- DimPlot(alpha, reduction ="pca", cols = c("darkblue"))
pca_p2 <- DimPlot(lambda, reduction ="pca", cols = c("blue"))
pca_p3 <- DimPlot(untreated, reduction ="pca", cols = c("lightblue"))

# Alternatively, the plots can be done together: 
grid.arrange(pca_p1, pca_p2, pca_p3, ncol = 3)

# tSNE 
alpha <- RunTSNE(alpha)
lambda <- RunTSNE(lambda)
untreated <- RunTSNE(untreated)

tsne_p1 <- DimPlot(alpha, reduction = "tsne", cols = "darkblue")
tsne_p2 <- DimPlot(lambda, reduction = "tsne", cols = "blue")
tsne_p3 <- DimPlot(untreated, reduction = "tsne", cols = "lightblue")
library(gridExtra)
grid.arrange(tsne_p1, tsne_p2, tsne_p3, ncol = 3)

# UMAP 
alpha <- RunUMAP(alpha, dims = 1:10)
lambda <- RunUMAP(lambda, dims = 1:10)
untreated <- RunUMAP(untreated, dims = 1:10)

umap_p1 <- DimPlot(alpha, reduction = "umap", cols = "darkblue")
umap_p2 <- DimPlot(lambda, reduction = "umap", cols = "blue")
umap_p3 <- DimPlot(untreated, reduction = "umap", cols = "lightblue")
library(gridExtra)
grid.arrange(umap_p1, umap_p2, umap_p3, ncol = 3)

# At this point I want to save the seurat ojects : 

##### Determine dimensionality of the dataset ####


alpha <- JackStraw(alpha, num.replicate = 100)
alpha <- ScoreJackStraw(alpha, dims = 1:20) # looks at scores for the first 20 PC
js_p1 <- JackStrawPlot(alpha, dims = 1:20)# Max dimension is 20


lambda <- JackStraw(lambda, num.replicate = 100)
lambda <- ScoreJackStraw(lambda, dims = 1:20)
js_p2 <- JackStrawPlot(lambda, dims = 1:20)

untreated <- JackStraw(untreated, num.replicate = 100)
untreated <- ScoreJackStraw(untreated, dims = 1:20)
js_p3 <- JackStrawPlot(untreated, dims = 1:20)

grid.arrange(js_p1, js_p2, js_p3, ncol = 3)

##### Cluster cells ##### 
# the first 20 principal components are being used to calculate the nearest neighbors for the cells in the alpha dataset.
library(ggplot2)
install.packages("RColorBrewer")
library(RColorBrewer)
palette.a <- brewer.pal(12, "Paired")
palette.b <- c(palette.a, "purple")

alpha <- FindNeighbors(alpha, dims = 1:20) 
lambda <- FindNeighbors(lambda, dims = 1:20)
untreated <- FindNeighbors(untreated, dims = 1:20)


alpha <- FindClusters(alpha, resolution = 0.5) 
lambda <- FindClusters(lambda, resolution = 0.5) 
untreated <- FindClusters(untreated, resolution = 0.5) 


dim_p1 <- DimPlot(alpha, group.by = "seurat_clusters", cols = palette.a) + ggtitle("IFN alpha")
dim_p2 <- DimPlot(lambda, group.by = "seurat_clusters", cols = palette.a)  + ggtitle("IFN lambda")
dim_p3 <- DimPlot(untreated, group.by = "seurat_clusters", cols = palette.b)  + ggtitle("Untreated")

grid.arrange(dim_p1, dim_p2, dim_p3, ncol = 3)

##### Finding differentially expressed features #####

alpha.markers <- FindMarkers(alpha, ident.1 = 1, min.pct = 0.25) # identifying markers and their differential expression information & store in seurat object
alpha <- AddMetaData(findmarkers = alpha.markers) # create a new slot named "markers" within the alpha Seurat object and assign the markers_result to it



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
is(alpha, "Seurat")
# [1] TRUE


# To launch cell type identification, we simply call the `classify_cells`function : 
alpha.scannotatR <- classify_cells(classify_obj = alpha, 
                             assay = 'RNA', 
                             slot = 'counts',
                             cell_types = 'all', 
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









##### Annotate clsters using SingleR #####
BiocManager::install("celldex")
BiocManager::install("ensembldb")

library(celldex)
library(SingleR)

ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl=TRUE) # reference dataset appropriate for PBMCs

# alpha.singleR <- alpha@assays$RNA@counts 

alpha.singleR <- as.SingleCellExperiment(alpha)
ref.singleR <- as.SingleCellExperiment(ref.data)

  
dim(alpha.singleR)
# [1] 17692  4742


dim(ref.data)
# 17893   713

alpha.pred <- SingleR(test = alpha.singleR,
                      ref = ref.data,
                      labels = ref.data$label.main)


# testing to see if making the objects have equal # of genes will fix it : 
# Assuming that 'common_genes' is a vector of gene names that are present in both datasets
common_genes <- intersect(rownames(alpha.singleR), rownames(ref.data))

# Subset both datasets to include only the common genes
alpha.singleR_subset <- alpha.singleR[common_genes, ]
ref.data_subset <- ref.data[common_genes, ]
dim(alpha.singleR_subset)
dim(ref.data_subset)

# fix : seeing if the inconsistent ggene naming could be the problem 
genes <- read_tsv("honours_2023/kallisto_index/transcripts_to_genes")
genes <- subset(genes, select = -ENST00000624431.2) # remove column 
genes$ENSG00000279928.2 <- substring(genes$ENSG00000279928.2, 1, 15) # take out version numbers 
colnames(genes) <- c("Ensembl_ID", "HGNC_ID") # changing column names 
genes$HGNC_ID <- ifelse(is.na(genes$HGNC_ID),genes$Ensembl_ID, genes$HGNC_ID)



# Now you can use the SingleR function with the subsetted datasets
alpha.pred <- SingleR(
  test = alpha.singleR_subset,
  ref = ref.data_subset,
  labels = ref.data$label.main
)



# output is cell barcodes as rows with columns for predictions 

# Now to visualize this, save the labels assigned by SingleR into the Seurat object 

alpha.pbmc.counts$singleR.labels <- alpha.pred$labels[match(rownames(alpha@meta.data), rownames(alpha.pred))]

# Plot this with labels : 

palette.a <- brewer.pal(11, "Paired")

Seurat::DimPlot(
  object = alpha,
  reduction = 'umap',
  group.by = 'singleR.labels',
  pt.size = 1,
  cols = palette.a
)


