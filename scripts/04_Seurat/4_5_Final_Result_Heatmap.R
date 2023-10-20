# Load required libraries
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")

library(reshape2)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)

TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated)

# [1] Find the average expression for a gene across treatment type: #####

Bcells <- subset(TreatmentAnnotated, celltype == "B")
result <- AverageExpression(Bcells, 
                  assay = "RNA",
                  features = c("CD79A","CD79B", "CD40", "IFNLR1","IFNAR1", # Receptors
                               "OAS1", "OAS2", "OAS3","IFIT1", "IFIT2", "IFIH1", "STAT1", # Core ISGs
                               "PAX5", "IRF9","IRF7",
                               "CD38", "JCHAIN", "XBP1","FAS", "CXCL10","BLNK", "CD19" # Differentiation markers 
                              # "", #Inflammation markers
                               ),
                  group.by = "treatment",
                  slot = "data",
                  verbose = TRUE)



result <- as.data.frame(result)
#result$gene <- rownames(result)
colnames(result)[1] <- "a"
colnames(result)[2] <- "l"
colnames(result)[3] <- "u"


# Convert the data frame to a matrix
result_matrix <- as.matrix(result)
col_names <- colnames(result)
row_names <- rownames(result)


# [2]  Heatmap #####

# COLOURS : "#95d16d", #2fa390, #5ac0d9, #719afb, #9a79c0

pheatmap(
  mat = result_matrix,
  labels_row = row_names,  # Row labels
  labels_col = col_names,
  cluster_rows = FALSE,    # Disable row clustering
  cluster_cols = TRUE,     # Enable column clustering
  color = colorRampPalette(c("#f1f3f3", "#95d16d"))(1000),  # Color palette
  treeheight_col = 0,  # Disable the dendrogram lines at the top
  display_numbers = FALSE   # Disable the display of numbers on the heatmap
)

# Define color palette
colors <- colorRampPalette(c("#f1f3f3", "#95d16d"))(100)

# Determine breaks and colors for the heatmap
heatmap_breaks <- seq(min(result_matrix), max(result_matrix), length.out = length(colors) + 1)
heatmap_colors <- cut(result_matrix, breaks = heatmap_breaks, labels = colors)

# Create circular heatmap
par(mar = c(1, 1, 1, 1))  # Set margins to fit the circular plot properly
circlize::circos.heatmap(result_matrix, col = colours, row_names = row_names, col_names = col_names)


# First seeing which genes are of interest : ######
# Feature Plots of individual gene expression 
FeaturePlot(TreatmentAnnotated, 
            features = c("CRP","SAA", "ICAM-1", "VCAM-1", "NFKB", "AP-1", "COX-2", "iNOS"),
            split.by = "sample",
            cols = c("grey", "black"))

Bcells <- subset(TreatmentAnnotated, celltype == "B")
DefaultAssay(Bcells) <- "RNA"

ISG20 <- FeaturePlot(Bcells, 
                     features = "ISG20",
                     split.by = "sample",
                     cols = c("grey", "black"),
                     pt.size = 5)

IFIT3 <- FeaturePlot(Bcells, 
                     features = "CD79B",
                     split.by = "treatment",
                     cols = c("lightgrey", "black"),
                     pt.size = 5) 

FeaturePlot(Bcells, 
            features = c("CD79A", "CD79B","CD40"),
            split.by = "sample",
            cols = c("grey", "black"),
            pt.size = 5) #+ 
# theme_classic()

ISG20 / IFIT3 / STAT1


CD38 <- FeaturePlot(Bcells, 
                    features = "CD38",
                    split.by = "",
                    cols = c("grey", "black"),
                    pt.size = 5)

JCHAIN <- FeaturePlot(Bcells, 
                      features = "JCHAIN",
                      split.by = "sample",
                      cols = c("grey", "black"),
                      pt.size = 5)

XBP1 <- FeaturePlot(Bcells, 
                    features = "XBP1",
                    split.by = "sample",
                    cols = c("grey", "black"),
                    pt.size = 5) #+ 
# theme_classic()

CD38 / JCHAIN / XBP1



?FeaturePlot
DefaultAssay(Bcells) <- "RNA"
Showme <- FeaturePlot(Bcells, 
                      features = c("CD79A", "CD79B", "CD40"), # , "CD79B", "IFNAR1", "IFN"
                      split.by = "treatment",
                      cols = c("grey", "black"),
                      pt.size = 5) 
Showme


DotPlot(Bcells, 
        features = c("mTORC1"), # CD79A", "CD79B "IFNLR1", "IL10RB", "TYK2""
        split.by = "sample",
        cols = c("black", "black", "black")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold")) + 
  coord_flip()


receptorDP <- DotPlot(Bcells,
                      features = c("CD40", "CXCR5"), #"CD79A", "CD79B", "IFNAR1", 
                      split.by = "sample",
                      cols = c("black", "black", "black")) + ##a9aaa9
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_text(face = "bold"), # angle = 90, hjust = 1, size = 14, 
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom") +
  geom_hline(yintercept = c(4.75, 5.25, 14.75, 15.25), color = "black", linetype = "dashed", size = 0.25) +
  coord_flip() +
  scale_y_discrete(labels = c("α", "λ", "u"))

receptorAllCells / receptorDP 
#| Showme





