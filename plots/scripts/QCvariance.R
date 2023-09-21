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


##### Seurat (dupe plot) ALPHA #####
# BEFORE 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(alpha@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(alpha@assays$RNA@counts, 1, var)
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
Log10Magnitude  <- log10(rowSums(alpha@assays$RNA@counts))

# calculate log10 of variance 
alpha <- FindVariableFeatures(alpha)
Log10Variance <- log10(alpha@assays$RNA@meta.features$vst.variance)

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

##### Seurat LAMBDA #####
# BEFORE 
# log10 magnitude for each cell
Log10Magnitude  <- log10(rowSums(lambda@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(lambda@assays$RNA@counts, 1, var)
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
Log10Magnitude  <- log10(rowSums(lambda@assays$RNA@counts))

# calculate log10 of variance 
#lambda <- FindVariableFeatures(lambda)
Log10Variance <- log10(lambda@assays$RNA@meta.features$vst.variance)

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
Log10Magnitude  <- log10(rowSums(untreated@assays$RNA@counts))

# calculate log10 of variance 
variances <- apply(untreated@assays$RNA@counts, 1, var)
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
Log10Magnitude  <- log10(rowSums(lambda@assays$RNA@counts))

# calculate log10 of variance 
#lambda <- FindVariableFeatures(lambda)
Log10Variance <- log10(lambda@assays$RNA@meta.features$vst.variance)

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
