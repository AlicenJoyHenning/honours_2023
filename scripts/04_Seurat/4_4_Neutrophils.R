# HEATMAP FOR VITAMIN D GENES IN NEUTROPHILS 

library(ggplot2)
library(pheatmap)
library(tibble)

TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated)
Neutrophils <- subset(TreatmentAnnotated, subset = celltype == "neutrophils") 

genes <- c(
           "CAMP", 
           "TSPAN18",
           "CD14", 
           "STEAP4",
           "CXCL10",   
           "S100A4", "S100A6","S100A8", "S100A9","S100A10", "S100A11", "S100A12", "S100A13", 
           "SRGN",
           "CD93", 
           "CEBPB", 
           "FN1", 
           "NINJ1", 
           "LILRB4", 
           "SEMA6B", 
           "THBD"
           )

avgexp <- AverageExpression(Neutrophils, assay = "RNA", return.seurat = T, group.by = 'treatment', colnames = TRUE)
logexp <- LogNormalize(Neutrophils)
avgexp <- AverageExpression(logexp, assay = "RNA", return.seurat = T, group.by = 'treatment', colnames = TRUE)


matrix <- avgexp@assays$RNA@counts
expression_df <- as.data.frame(matrix, rownames = TRUE)
expression_df <- rownames_to_column(expression_df, var = "Gene")
expression_df <- subset(expression_df, Gene %in% genes)

matrix_data <- as.matrix(expression_df[, -1])
row_names <- expression_df$Gene        


pheatmap(
  mat = matrix_data,
  labels_row = row_names,  # Row labels
  cluster_rows = FALSE,    # Disable row clustering
  cluster_cols = TRUE,     # Enable column clustering
  color = colorRampPalette(c("#6ab5ba", "white", "grey"))(50),  # Color palette
  main = "Average Gene Expression Heatmap"
)

##### Seurat Default heatmap #####
DoHeatmap(
  object = Neutrophils,
  features = genes,
 # cells = NULL,
  group.by = "treatment",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = NULL,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
)
  
  
