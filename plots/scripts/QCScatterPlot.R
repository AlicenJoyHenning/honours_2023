# Scatter plots for cell quality 
install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(ggpubr)


##### ALPHA #####
Amatrix <- "honours/work/1109/alpha/matrix.mtx.gz"
Abarcodes <- "honours/work/1109/alpha/barcodes.tsv.gz"
Afeatures <- "honours/work/1109/alpha/features.tsv.gz"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=0) 


alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")


# # Assuming that alpha is a Seurat object and nFeature_RNA, nCount_RNA, and quality are columns within the meta.data slot
alpha@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA >= 2500 & alpha@meta.data$percent.mt < 10] <- "HIGH nF"
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA < 2500] <- "Accepted"
alpha@meta.data$Quality[alpha@meta.data$nFeature_RNA >= 2500 & alpha@meta.data$percent.mt >= 10] <- "Accepted"

colours <- c("lightgrey","#696969")

# Scatter for features 
alphascatter <- 
  ggplot(alpha@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
   # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "α",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7))) # Adjust point size in the legend

# Scatter for mitochondrial percent 

alpha@meta.data$MQuality <- NA  # Initialize a quality column in the meta.data slot

alpha@meta.data$MQuality[alpha@meta.data$percent.mt <= 10] <- "Accepted"
alpha@meta.data$MQuality[alpha@meta.data$percent.mt > 10] <- "High Mt"

mcolours <- c("lightgrey", "#6ab5ba")

alphamscatter <- 
  ggplot(alpha@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = MQuality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = mcolours) + 
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "α",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7)))

##### LAMBDA #####
Lmatrix <- "honours/work/1109/lambda/matrix.mtx.gz"
Lbarcodes <- "honours/work/1109/lambda/barcodes.tsv.gz"
Lfeatures <- "honours/work/1109/lambda/features.tsv.gz"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features = 0)

# quality columns within the meta.data slot
lambda@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA >= 2500 & lambda@meta.data$percent.mt < 10] <- "HIGH nF"
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA < 2500] <- "Accepted"
lambda@meta.data$Quality[lambda@meta.data$nFeature_RNA >= 2500 & lambda@meta.data$percent.mt >= 10] <- "Accepted"

colours <- c("lightgrey","#696969")


# Scatter for features 
lambdascatter <- 
  ggplot(lambda@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "λ",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7))) # Adjust point size in the legend

# Scatter for mitochondrial percent 
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
lambda@meta.data$MQuality <- NA  # Initialize a quality column in the meta.data slot

lambda@meta.data$MQuality[lambda@meta.data$percent.mt <= 10] <- "Accepted"
lambda@meta.data$MQuality[lambda@meta.data$percent.mt > 10] <- "High Mt"

mcolours <- c("lightgrey", "#6ab5ba")

lambdamscatter <- 
  ggplot(lambda@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = MQuality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = mcolours) + 
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_blank(),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "λ",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7)))


##### UNTREATED ####
Umatrix <- "honours/work/1109/untreated/matrix.mtx.gz"
Ubarcodes <- "honours/work/1109/untreated/barcodes.tsv.gz"
Ufeatures <- "honours/work/1109/untreated/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features = 0)

# quality columns within the meta.data slot
untreated@meta.data$Quality <- NA  # Initialize a quality column in the meta.data slot

# Assign quality values based on nFeature_RNA values
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA >= 2500 & untreated@meta.data$percent.mt < 10] <- "HIGH nF"
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA < 2500] <- "Accepted"
untreated@meta.data$Quality[untreated@meta.data$nFeature_RNA >= 2500 & untreated@meta.data$percent.mt >= 10] <- "Accepted"

colours <- c("lightgrey","#696969")

# Scatter for features 
untreatedscatter <- 
  ggplot(untreated@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = Quality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = colours) + 
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_text(size = 16, face = "bold", vjust = 0.01),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "u",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7))) # Adjust point size in the legend

# Scatter for mitochondrial percent 
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")
untreated@meta.data$MQuality <- NA  # Initialize a quality column in the meta.data slot

untreated@meta.data$MQuality[untreated@meta.data$percent.mt <= 10] <- "Accepted"
untreated@meta.data$MQuality[untreated@meta.data$percent.mt > 10] <- "High Mt"

mcolours <- c("lightgrey", "#6ab5ba")

untreatedmscatter <- 
  ggplot(untreated@meta.data, aes(x = nFeature_RNA, y = nCount_RNA, color = MQuality)) +
  geom_point(size = 1) +  # Set the size to a constant value (e.g., 3)
  scale_color_manual(values = mcolours) + 
  theme_classic() +
  theme(
    # panel.background = element_rect(fill = "#d3d3d3"),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16, face = "bold", vjust = 0.01),
    axis.title.y = element_text(size = 16, face = "bold", vjust = 0.01),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  labs(title = "u",
       x = "Features",
       y = "UMIs") +
  guides(color = guide_legend(override.aes = list(size = 7)))

##### combine #####


(untreatedscatter | lambdascatter | alphascatter ) / (untreatedmscatter | lambdamscatter | alphamscatter)


