# MANUAL WILCOXON RANK SUM TEST FOR GENES NOT IDENTIFIED BY FINDCONSERVEDMARKERS

##### Load dataset and create testing subsets for each cluster ####
treatment <- readRDS("honours/work/1109/treatment.rds")

# For each cluster, we create 2 groups: 1 containing cells belonging to the
# cluster of interest, and the other containing cells that do not belong to the
# cluster of interest (aka all the other clusters)

treatment@meta.data$samples <- treatment@meta.data$treatment

##### Run individual Wilcoxon Rank sum tests for genes of interest ####
##### CLUSTER 9 NKT cells ####
cluster9 <- subset(treatment, subset = seurat_clusters == 9)
notcluster9 <- subset(treatment, subset = seurat_clusters != 9)

# CD8A
CD8ACluster9 <- t(cluster9@assays$integrated@data["CD8A", ])
CD8ACluster9 <- data.frame(Cell = colnames(CD8ACluster9),Expression = CD8ACluster9[,])
CD8AnotCluster9 <- t(notcluster9@assays$integrated@data["CD8A", ])
CD8AnotCluster9 <- data.frame(Cell = colnames(CD8AnotCluster9),Expression = CD8AnotCluster9[,])
CD8Awilcox_result <- wilcox.test(CD8ACluster9$Expression, CD8AnotCluster9$Expression)
dim <- dim(cluster11@assays$integrated@data) # find # comparisons 
dim <- dim[1]
CD8Aresult <- CD8Awilcox_result$p.value * dim
print(CD8Aresult)

# NK1.1 : not present in data set
NK1Cluster9 <- t(cluster9@assays$integrated@data["NK1.1", ])
Nk1Cluster9 <- data.frame(Cell = colnames(NK1Cluster9),Expression = NK1Cluster9[,])
NK1notCluster9 <- t(notcluster9@assays$integrated@data["NK1.1", ])
NK1AnotCluster9 <- data.frame(Cell = colnames(NK1notCluster9),Expression = NK1AnotCluster9[,])
NK1wilcox_result <- wilcox.test(NK1Cluster9$Expression, NK1AnotCluster9$Expression)
dim <- dim(cluster9@assays$integrated@data) # find # comparisons 
dim <- dim[1]
NK1result <- NK1wilcox_result$p.value * dim
print(NK1result)

# CD8B
CD8BCluster9 <- t(cluster9@assays$integrated@data["CD8B", ])
CD8BCluster9 <- data.frame(Cell = colnames(CD8BCluster9),Expression = CD8BCluster9[,])
CD8BnotCluster9 <- t(notcluster9@assays$integrated@data["CD8B", ])
CD8BnotCluster9 <- data.frame(Cell = colnames(CD8BnotCluster9),Expression = CD8BnotCluster9[,])
CD8Bwilcox_result <- wilcox.test(CD8BCluster9$Expression, CD8BnotCluster9$Expression)
dim <- dim(cluster9@assays$integrated@data) # find # comparisons 
dim <- dim[1]
CD8Bresult <- CD8Bwilcox_result$p.value * dim
print(CD8Bresult)

# GZMA
GZMACluster9 <- t(cluster9@assays$integrated@data["BLNK", ])
GZMACluster9 <- data.frame(Cell = colnames(GZMACluster9),Expression = GZMACluster9[,])
GZMAnotCluster9 <- t(notcluster9@assays$integrated@data["BLNK", ])
GZMAnotCluster9 <- data.frame(Cell = colnames(GZMAnotCluster9),Expression = GZMAnotCluster9[,])
GZMAwilcox_result <- wilcox.test(GZMACluster9$Expression, GZMAnotCluster9$Expression)
dim <- dim(cluster9@assays$integrated@data) # find # comparisons 
dim <- dim[1]
GZMAresult <- GZMAwilcox_result$p.value * dim
print(GZMAresult)

##### CLUSTER 10 T regs #####

cluster10 <- subset(treatment, subset = seurat_clusters == 10)
notcluster10 <- subset(treatment, subset = seurat_clusters != 10) 
# can test in individual treatment conditions using the following : 
notcluster10 <- subset(treatment, subset = seurat_clusters != 10 & samples == "alpha") 

# FOXP3 
# Extract the expression values for the gene of interest in the two groups
FOXP3Cluster10 <- t(cluster10@assays$RNA@data["FOXP3", ])
FOXP3Cluster10 <- data.frame(
  Cell = colnames(FOXP3Cluster10),
  Expression = FOXP3Cluster10[,]
  )

FOXP3notCluster10 <- t(notcluster10@assays$RNA@data["FOXP3", ])
FOXP3notCluster10 <- data.frame(
  Cell = colnames(FOXP3notCluster10),
  Expression = FOXP3notCluster10[,]
)

# Perform the Wilcoxon Rank Sum test
FOXP3wilcox_result <- wilcox.test(FOXP3Cluster10$Expression,FOXP3notCluster10$Expression)

dim <- dim(notcluster10@assays$RNA@data) # find # comparisons 
dim <- dim[1]
# bonferroni correction : 
FOXP3result <- FOXP3wilcox_result$p.value * dim
print(FOXP3result)

# CCR4 
CCR4Cluster10 <- t(cluster10@assays$RNA@data["CCR4", ])
CCR4Cluster10 <- data.frame(Cell = colnames(CCR4Cluster10),Expression = CCR4Cluster10[,])
CCR4notCluster10 <- t(notcluster10@assays$RNA@data["CCR4", ])
CCR4notCluster10 <- data.frame(Cell = colnames(CCR4notCluster10),Expression = CCR4notCluster10[,])
CCR4wilcox_result <- wilcox.test(CCR4Cluster10$Expression, CCR4notCluster10$Expression)
dim <- dim(notcluster10@assays$RNA@data) # find # comparisons 
dim <- dim[1]
CCR4result <- CCR4wilcox_result$p.value * dim
print(CCR4result)

# IKZF2 
IKZF2Cluster10 <- t(cluster10@assays$RNA@data["IKZF2", ])
IKZF2Cluster10 <- data.frame(Cell = colnames(IKZF2Cluster10), Expression = IKZF2Cluster10[,])
IKZF2notCluster10 <- t(notcluster10@assays$RNA@data["IKZF2", ])
IKZF2notCluster10 <- data.frame(Cell = colnames(IKZF2notCluster10),Expression = IKZF2notCluster10[,])
IKZF2wilcox_result <- wilcox.test(IKZF2Cluster10$Expression, IKZF2notCluster10$Expression)
dim <- dim(cluster10@assays$RNA@data) # find # comparisons 
dim <- dim[1]
IKZF2result <- (IKZF2wilcox_result$p.value) * (dim)
print(IKZF2result)



# FANK1
FANK1Cluster10 <- t(cluster10@assays$RNA@data["FANK1", ])
FANK1Cluster10 <- data.frame(Cell = colnames(FANK1Cluster10), Expression = FANK1Cluster10[,])
FANK1notCluster10 <- t(notcluster10@assays$RNA@data["FANK1", ])
FANK1notCluster10 <- data.frame(Cell = colnames(FANK1notCluster10), Expression = FANK1notCluster10[,])
FANK1wilcox_result <- wilcox.test(FANK1Cluster10$Expression, FANK1notCluster10$Expression)
dim <- dim(cluster10@assays$RNA@data) # find # comparisons 
dim <- dim[1]
FANK1result <- FANK1wilcox_result$p.value * dim
print(FANK1result)


##### CLUSTER 11 NK cells ####
cluster11 <- subset(treatment, subset = seurat_clusters == 11)
notcluster11 <- subset(treatment, subset = seurat_clusters != 11) 

# GNLY 
GNLYCluster11 <- t(cluster11@assays$integrated@data["GNLY", ])
GNLYCluster11 <- data.frame(Cell = colnames(GNLYCluster11),Expression = GNLYCluster11[,])
GNLYnotCluster11 <- t(notcluster11@assays$integrated@data["GNLY", ])
GNLYnotCluster11 <- data.frame(Cell = colnames(GNLYnotCluster11),Expression = GNLYnotCluster11[,])
GNLYwilcox_result <- wilcox.test(GNLYCluster11$Expression, GNLYnotCluster11$Expression)
dim <- dim(cluster11@assays$integrated@data) # find # comparisons 
dim <- dim[1]
GNLYresult <- GNLYwilcox_result$p.value * dim
print(GNLYresult)
GNLYwilcox_result$p.value

# GZMB 
GZMBCluster11 <- t(cluster11@assays$integrated@data["GZMB", ])
GZMBCluster11 <- data.frame(Cell = colnames(GZMBCluster11),Expression = GZMBCluster11[,])
GZMBnotCluster11 <- t(notcluster11@assays$integrated@data["GZMB", ])
GZMBnotCluster11 <- data.frame(Cell = colnames(GZMBnotCluster11),Expression = GZMBnotCluster11[,])
GZMBwilcox_result <- wilcox.test(GZMBCluster11$Expression, GZMBnotCluster11$Expression)
dim <- dim(cluster11@assays$integrated@data) # find # comparisons 
dim <- dim[1]
GZMBresult <- GZMBwilcox_result$p.value * dim
print(GZMBresult)
GZMBwilcox_result$p.value

# GZMBCluster11 <- t(cluster11@assays$RNA@data["GZMB", ])
NKG7Cluster11 <- t(cluster11@assays$integrated@data["NKG7", ])
NKG7Cluster11 <- data.frame(Cell = colnames(NKG7Cluster11),Expression = NKG7Cluster11[,])
NKG7notCluster11 <- t(notcluster11@assays$integrated@data["NKG7", ])
NKG7notCluster11 <- data.frame(Cell = colnames(NKG7notCluster11),Expression = NKG7notCluster11[,])
NKG7wilcox_result <- wilcox.test(NKG7Cluster11$Expression, NKG7notCluster11$Expression)
dim <- dim(cluster11@assays$integrated@data) # find # comparisons 
dim <- dim[1]
NKG7result <- NKG7wilcox_result$p.value * dim
print(NKG7result)


