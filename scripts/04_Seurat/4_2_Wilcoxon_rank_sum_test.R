# MANUAL WILCOXON RANK SUM TEST FOR GENES NOT IDENTIFIED BY FINDCONSERVEDMARKERS

##### [1] Load dataset and create testing subsets for each cluster ####
treatment <- readRDS("honours/work/1109/treatment.rds")
DefaultAssay(treatment) <- "RNA"
colnames(treatment@meta.data)[1] <- "sample"


# For each cluster, we create 2 groups: 1 containing cells belonging to the
# cluster of interest, and the other containing cells that do not belong to the
# cluster of interest (aka all the other clusters)

##### [2] Function #####

  
conditions <- {c("alpha", "lambda", "untreated")
genes0 <- c("CD14", "FCGR3A", "CSF3R", "CXCL8", "ITGAX") # cluster 0 and 4
genes1 <- c("LEF1", "CD28") # cluster 1
genes2 <- c("CD14", "VASP", "SOD2", "SPI1", "JAML") # cluster 2 
genes3 <- c("IL7R", "CD2","IL32", "GPR183") # T cells
genes5 <- c("CD8A", "CD8B", "IL7R", "GAS5","PASK") # naive
genes6 <- c("CD37", "CD24", "BCL11A", "BLNK") # B 
genes7 <- c("CD8A", "CD8B", "CD3G", "GZMA", "GZMB") # NKT
genes8 <- c("S100A12", "LMNA","CD14", "MS4A7", "ITGAX", "FUT4") # mDCs
genes9 <- c("CD8A", "CD8B","GZMA", "CD3G") # cytotoxic T cells 
genes10 <- c("FOXP3", "CCR4","FANK1") # T regs 
genes11 <- c("GNLY", "GZMB","NKG7", "CD3G") # NK
genes12 <- c("CLEC1B", "GP9","PF4", "PPBP") # platelets
genes14 <- c("ITGAX", "CCRL2", "CTSD", "GP9") # DCs
genes15 <- c("CD3G", "NRGN", "CCR7") # CD4 T 
genes17 <- c("CLEC4C", "GZMB","GAPT", "IRF7", "PLD4") # pDCs
}

DIYWilcoxon <- function(condition, CellType, gene) {
  print(paste("Performing Wilcoxon Rank Sum test for gene", gene, "in ", CellType))
  
  # Create subset for the specified condition
  conditionSubset <- subset(TreatmentAnnotated, subset = sample == condition)
  
  # Create subsets for the specified cluster and its complement
  clusterSubset <- subset(conditionSubset, subset = celltype == CellType)
  notClusterSubset <- subset(conditionSubset, subset = celltype != CellType)
  
  # Extract the expression values for the specified gene in the two groups
  geneCluster <- t(clusterSubset@assays$RNA@data[gene, ])
  geneCluster <- data.frame(
    Cell = colnames(geneCluster),
    Expression = geneCluster[,]
  )
  
  geneNotCluster <- t(notClusterSubset@assays$RNA@data[gene, ])
  geneNotCluster <- data.frame(
    Cell = colnames(geneNotCluster),
    Expression = geneNotCluster[,]
  )
  
  # Perform the Wilcoxon Rank Sum test
  wilcoxResult <- wilcox.test(geneCluster$Expression, geneNotCluster$Expression)
  
  # Bonferroni correction
  dim <- dim(notClusterSubset@assays$integrated@data)[1]  # find # comparisons 
  adjResult <- wilcoxResult$p.value * dim
  
  print(paste("Adj P value", "|", condition, "|", adjResult))
  return(adjResult)
}

for (gene in genes1) {
  DIYWilcoxon("alpha", "", gene)
}











##### [*] Preliminary testing for individual Wilcoxon Rank sum calculations for genes of interest ####
# Remove redundant phrasing with 'treatment' (Seurat object & meta data column)
colnames(treatment@meta.data)[1] <- "sample"

# Create subset for just alpha 
alphasamples <- subset(treatment, subset = sample == "alpha")

# Create subset for cluster 10 
cluster10 <- subset(alphasamples, subset = seurat_clusters == 10 )
notcluster10 <- subset(alphasamples, subset = seurat_clusters != 10) 

# FOXP3 
# Extract the expression values for the gene of interest in the two groups
FOXP3Cluster10 <- t(cluster10@assays$RNA@data["BANK1", ])
FOXP3Cluster10 <- data.frame(
  Cell = colnames(FOXP3Cluster10),
  Expression = FOXP3Cluster10[,]
)

FOXP3notCluster10 <- t(notcluster10@assays$RNA@data["BANK1", ])
FOXP3notCluster10 <- data.frame(
  Cell = colnames(FOXP3notCluster10),
  Expression = FOXP3notCluster10[,]
)

# Perform the Wilcoxon Rank Sum test
FOXP3wilcox_result <- wilcox.test(FOXP3Cluster10$Expression,FOXP3notCluster10$Expression)

dim <- dim(notcluster10@assays$integrated@data) # find # comparisons 
dim <- dim[1]
# bonferroni correction : 
FOXP3result <- FOXP3wilcox_result$p.value * dim
print(FOXP3result)


