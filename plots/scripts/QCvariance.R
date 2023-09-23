# Scatter plots for normalization and variance visualisation 

BiocManager::install("pagoda2")

library(Matrix)
library(dplyr)
library(ggplot2)
library(igraph)
library(pagoda2)
library(Seurat)
library(SeuratObject)

##### Pagoda2 adjustVariance() plot #####

path <- "honours/work/1109/alpha/"
alpha <- read10xMatrix(path, version = "V3", transcript.id = "SYMBOL", verbose = TRUE)
dim(alpha)

counts <- gene.vs.molecule.cell.filter(alpha, min.cell.size=200, max.cell.size = 2500)
rownames(counts) <- make.unique(rownames(counts))
alpha <- Pagoda2$new(counts, log.scale=TRUE, n.cores=1)

alpha$adjustVariance(plot=TRUE, gam.k=10)

# although the output is beautiful, it won't allow me to use the plotting function without performing their variance methods, I wanted to use it to visualize the effects from Seurat 


##### Perform filtering #####
# Load datasets: alpha, lambda, and untreated using ReadMtx function : 
Amatrix <- "honours/work/1109/alpha/matrix.mtx.gz"
Abarcodes <- "honours/work/1109/alpha/barcodes.tsv.gz"
Afeatures <- "honours/work/1109/alpha/features.tsv.gz"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=0)

Lmatrix <- "honours/work/1109/lambda/matrix.mtx.gz"
Lbarcodes <- "honours/work/1109/lambda/barcodes.tsv.gz"
Lfeatures <- "honours/work/1109/lambda/features.tsv.gz"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features = 0)

Umatrix <- "honours/work/1109/untreated/matrix.mtx.gz"
Ubarcodes <- "honours/work/1109/untreated/barcodes.tsv.gz"
Ufeatures <- "honours/work/1109/untreated/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features = 0)


alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

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

##### Seurat (dupe plot) ALPHA #####
# BEFORE 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(alpha@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(alpha@assays$RNA@counts, 1, var)
Log10Variance <- log10(variances)

df <- data.frame(Log10Magnitude, Log10Variance)

before <- 
  ggplot(data = df, 
       aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "before") +
  guides(color = guide_legend(override.aes = list(size = 7)))


# AFTER 
Log10Magnitude  <- log10(rowSums(alpha@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(alpha@assays$RNA@counts, 1, var)
Log10Variance <- log10(variances)

df <- data.frame(Log10Magnitude, Log10Variance)


# log10 magnitude for each cell
# Log10Magnitude  <- log10(rowSums(alphaAF@assays$RNA@counts))
# # calculate log10 of variance 
# alpha <- FindVariableFeatures(alphaAF)
# Log10Variance <- log10(alphaAF@assays$RNA@meta.features$vst.variance)

# plot 
after <- 
  ggplot(data = df, 
         aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "adjusted") +
  guides(color = guide_legend(override.aes = list(size = 7)))

before | after 

##### Seurat LAMBDA #####
# BEFORE 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(lambdaNF@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(lambdaNF@assays$RNA@counts, 1, var)
Log10Variance <- log10(variances)


df <- data.frame(Log10Magnitude, Log10Variance)

# plot 

before <- 
  ggplot(data = df, 
         aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "before") +
  guides(color = guide_legend(override.aes = list(size = 7)))


# AFTER 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(lambdaAF@assays$RNA@counts))

# calculate log10 of variance 
#lambda <- FindVariableFeatures(lambda)
Log10Variance <- log10(lambdaAF@assays$RNA@meta.features$vst.variance)

df <- data.frame(Log10Magnitude, Log10Variance)

# plot 
after <- 
  ggplot(data = df, 
         aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "adjusted") +
  guides(color = guide_legend(override.aes = list(size = 7)))

before | after 

##### Seurat UNTREATED ####
# BEFORE 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(untreatedBF@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(untreatedBF@assays$RNA@counts, 1, var)
Log10Variance <- log10(variances)


df <- data.frame(Log10Magnitude, Log10Variance)

# plot 

before <- 
  ggplot(data = df, 
         aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "before") +
  guides(color = guide_legend(override.aes = list(size = 7)))


# AFTER 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(lambdaAF@assays$RNA@counts))

# calculate log10 of variance 
#lambda <- FindVariableFeatures(lambda)
Log10Variance <- log10(lambdaAF@assays$RNA@meta.features$vst.variance)

df <- data.frame(Log10Magnitude, Log10Variance)

# plot 
after <- 
  ggplot(data = df, 
         aes(x = Log10Magnitude, y = Log10Variance)) +
  geom_point(colour = "#c9cacd") +
  labs(x = "log10(magnitude)", y = "log10(variance)") +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  labs(title = "adjusted") +
  guides(color = guide_legend(override.aes = list(size = 7)))

before | after 
