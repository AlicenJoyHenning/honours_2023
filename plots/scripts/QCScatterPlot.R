# Scatter plot for cell quality 
install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(ggpubr)


##### ALPHA #####

AlphaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/matrix.mtx.gz","honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)

alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")

# Assuming that alpha is a Seurat object and nFeature_RNA, nCount_RNA, and quality are columns within the meta.data slot
alpha@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA < 200] <- "LOW nF"
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA > 2500] <- "HIGH nF"
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA >= 200 & alpha@meta.data$nFeature_RNA <= 2500] <- "Accepted"

colours <- c("#6ab5ba", "white")

alphascatter <- 
  ggplot(alpha@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "#d3d3d3"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "alpha") +
  guides(color = guide_legend(override.aes = list(size = 7))) # Adjust point size in the legend
  

##### LAMBDA #####

LambdaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)

lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")

lambda@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA < 200] <- "LOW nF"
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA > 2500] <- "HIGH nF"
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA >= 200 & lambda@meta.data$nFeature_RNA <= 2500] <- "Accepted"

lambdascatter <- 
  ggplot(lambda@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "#d3d3d3"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "lambda") +
  guides(color = guide_legend(override.aes = list(size = 7))) 

##### UNTREATED #####
UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz")
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

untreated@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA < 200] <- "LOW nF"
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA > 2500] <- "HIGH nF"
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA >= 200 & untreated@meta.data$nFeature_RNA <= 2500] <- "Accepted"

untreatedscatter <- 
  ggplot(untreated@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "#d3d3d3"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "untreated") +
  guides(color = guide_legend(override.aes = list(size = 7))) 

##### combine #####

ggarrange(alphascatter, lambdascatter, untreatedscatter, ncol = 3)
remove_legend(alphascatter)

alphascatter | lambdascatter | untreatedscatter


