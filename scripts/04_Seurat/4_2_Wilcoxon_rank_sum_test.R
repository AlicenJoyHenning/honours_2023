# MANUAL WILCOXON RANK SUM TEST FOR GENES NOT IDENTIFIED BY FINDCONSERVEDMARKERS

##### [1] Load dataset and create testing subsets for each cluster ####
treatment <- readRDS("honours/work/1109/treatment.rds")
DefaultAssay(treatment) <- "RNA"

# For each cluster, we create 2 groups: 1 containing cells belonging to the
# cluster of interest, and the other containing cells that do not belong to the
# cluster of interest (aka all the other clusters)

##### [1.2] Test for individual Wilcoxon Rank sum tests for genes of interest ####
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



##### [2] Function #####
  
DIYWilcoxon <- function(condition, cluster, gene) {
  print(paste("Performing Wilcoxon Rank Sum test for gene", gene, "in Cluster", cluster, "for condition", condition, ":\n"))
  
  # Create subset for the specified condition
  conditionSubset <- subset(treatment, subset = sample == condition)
  
  # Create subsets for the specified cluster and its complement
  clusterSubset <- subset(conditionSubset, subset = seurat_clusters == cluster)
  notClusterSubset <- subset(conditionSubset, subset = seurat_clusters != cluster)
  
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
  
  print(paste("Adjusted P value for", gene, ":", adjResult, "\n"))
  return(adjResult)
}

DIYWilcoxon("alpha", 10, "BANK1")  


  













