# Sequencing Depth Boxplot 

library(Seurat) # using Seurat object 
library(Matrix)
library(ggplot2)

setwd("")
getwd()


# Load datasets: alpha, lambda, and untreated using ReadMtx function : 
Amatrix <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/unfiltered_counts.mtx"
Abarcodes <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/barcodes.tsv.gz"
Afeatures <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/features.tsv.gz"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=0)
saveRDS(alpha, "honours/work/1109/alpha/alphaCountMatrix.rds")

Lmatrix <- "honours/work/1109/lambda/matrix.mtx.gz"
Lbarcodes <- "honours/work/1109/lambda/barcodes.tsv.gz"
Lfeatures <- "honours/work/1109/lambda/features.tsv.gz"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features = 0)
saveRDS(lambda, "honours/work/1109/lambda/lambdaCountMatrix.rds")

Umatrix <- "honours/work/1109/untreated/matrix.mtx.gz"
Ubarcodes <- "honours/work/1109/untreated/barcodes.tsv.gz"
Ufeatures <- "honours/work/1109/untreated/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features = 0)
saveRDS(untreated, "honours/work/1109/untreated/untreatedCountMatrix.rds")

Aunfiltered <- Read10X("D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/")

?Read10X()
