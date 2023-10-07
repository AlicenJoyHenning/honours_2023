# HEATMAP FOR VITAMIN D GENES IN NEUTROPHILS 

library(ggplot2)
library(pheatmap)
library(tibble)

install.packages("psych") # finding geometric mean 
library(psych)

##### DIY #####
TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated)

# Subset to isolate neutrophils : 
Neutrophils <- subset(TreatmentAnnotated, subset = celltype == "neutrophils") 

# Create count matrices for each treatment type : 
NAlpha <- subset(Neutrophils, subset = treatment == "alpha") 
NMAlpha <- as.matrix(NAlpha@assays$RNA@counts)
NMAlpha <- log2(NMAlpha + 1) 
# Alpha <- apply(NMAlpha, 1, median)
Alpha <- apply(NMAlpha, 1, max) # Find the minimum value for each row 
Alpha <- as.matrix(Alpha)
dim(Alpha)

NLambda <- subset(Neutrophils, subset = treatment == "lambda") 
NMLambda <- as.matrix(NLambda@assays$RNA@counts)
NMLambda <- log2(NMLambda + 1)
# Lambda <- apply(NMLambda, 1, mean)
Lambda <- apply(NMLambda, 1, max)
Lambda <- as.matrix(Lambda)
dim(Lambda)


NUntreated <- subset(Neutrophils, subset = treatment == "untreated") 
NMUntreated <- as.matrix(NUntreated@assays$RNA@counts)
NMUntreated <- log2(NMUntreated + 1)
Untreated <- apply(NMUntreated, 1, max)
Untreated <- as.matrix(Untreated)
dim(Untreated)

# combine the median log values for each datasets for each gene into one matrix
NeutrophilMatrix <- cbind(Alpha, Lambda, Untreated)
dim(NeutrophilMatrix)

genes <- c("CAMP", "CD14", "CXCL10", "CD93", "LILRB4", "SEMA6B", "THBD", "IL10", "CCL2", "DEFA1")


Neutrophildf <- as.data.frame(NeutrophilMatrix)
Neutrophildf <- rownames_to_column(Neutrophildf, var = "Gene")
Neutrophildf <- subset(Neutrophildf, Gene %in% genes)
colnames(Neutrophildf) <- c("gene", "alpha", "lambda", "untreated")

FinalNeutrophilMatrix <- as.matrix(Neutrophildf[, -1]) # taking away gene column
row_names <- Neutrophildf$gene # this is where we keep the info in the gene column for labelling 


pheatmap(
  mat = FinalNeutrophilMatrix,
  labels_row = row_names,  # Row labels
  cluster_rows = FALSE,    # Disable row clustering
  cluster_cols = TRUE,     # Enable column clustering
  color = colorRampPalette(c("white", "#6ab5ba"))(5000),  # Color palette
  main = "Neutrophil Gene expression \n (log2(n + 1)) \n following IFN treatment"
)



### 
# genes <- c("CAMP", 
#            "TSPAN18",
#            "CD14", 
#            "STEAP4",
#            "CXCL10",   
#            "S100A4", "S100A6","S100A8", "S100A9","S100A10", "S100A11", "S100A12", "S100A13", 
#            "SRGN",
#            "CD93", 
#            "CEBPB", 
#            "FN1", 
#            "NINJ1", 
#            "LILRB4", 
#            "SEMA6B", 
#            "THBD"
#            )

genes <- c("CAMP", "CD14", "CXCL10", "CD93", "LILRB4", "SEMA6B", "THBD", "IL10", "CCL2", "DEFA1")

avgexp <- AverageExpression(Neutrophils, assay = "RNA", group.by = 'treatment', colnames = TRUE)
# logexp <- LogNormalize(Neutrophils)
# avgexp <- AverageExpression(logexp, assay = "RNA", return.seurat = T, group.by = 'treatment', colnames = TRUE)


avgexp <- as.data.frame(avgexp, rownames = TRUE)
avgexp <- rownames_to_column(avgexp, var = "Gene")
avgexp <- subset(avgexp, Gene %in% genes)

matrix_data <- as.matrix(avgexp[, -1])
row_names <- avgexp$Gene        

pheatmap(
  mat = matrix_data,
  labels_row = row_names,  # Row labels
  cluster_rows = FALSE,    # Disable row clustering
  cluster_cols = TRUE,     # Enable column clustering
  color = colorRampPalette(c("grey", "white", "#6ab5ba"))(50),  # Color palette
  annotation_col = NULL,   # Remove dendrogram lines at the top
  annotation_row = NULL,   # Remove dendrogram lines on the left
  main = "Average Gene Expression Heatmap",
  annotation_names_col = FALSE,  # Remove column annotation labels
  annotation_names_row = FALSE   # Remove row annotation labels
)







##### Seurat Default heatmap #####
genes <- c(
  "CXCL10",  
  "CD93",  
  "LILRB4", 
  "SEMA6B", 
  "THBD",
  "CD14",
  "CD38"
)

genes <- c(
  "TLR2",
  "NOD2",
  "CCL2",
  "IL10",
  "IL1RN",
  "CCL5",
  "CXCL10"
  )

genes <- c(
  "ITL3"
)
  


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
  
###### Gene sets #####
neutrophilDEGs <- readRDS("honours/results/FinalIndex/DGEAnalysis/RNAAssayDEGs/NeutrophilsLambdaResponse.rds")
genes <- neutrophilDEGs$gene
genes <- genes[-2]
