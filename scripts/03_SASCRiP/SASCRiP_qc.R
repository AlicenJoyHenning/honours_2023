# TESTING VALIDITY OF RESULTS FOLLOWING SASCRiP 

##### [1] Load dependencies #####
getwd()
setwd("~")

library(BiocManager)
BiocManager::install('limma') # BiocManager::install('limma')
BiocManager::install("SeuratData")
BiocManager::install('multtest')
BiocManager::install('metap')
BiocManager::install('xlsx')
BiocManager::install('pheatmap')
BiocManager::install('XLConnect')
BiocManager::install('writexl')


library(dplyr)
library(writexl)
library(openxlsx)
library(XLConnect)
library(ggplot2)
library(grid)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(SingleR)
library(celldex)
library(patchwork)
library(readr)
library(Matrix)
library(metap)
library(openxlsx)
library(cowplot)
library(readxl)
library(xlsx)

##### [2] Load datasets: alpha, lambda, and untreated using ReadMtx function : ####

Amatrix <- "../../media/intel6700/Passport/AlicenHonours/work/0909/alpha/seurat_matrix/matrix.mtx.gz"
Abarcodes <- "../../media/intel6700/Passport/AlicenHonours/work/0909/alpha/seurat_matrix/barcodes.tsv.gz"
Afeatures <- "../../media/intel6700/Passport/AlicenHonours/work/0909/alpha/seurat_matrix/features.tsv.gz"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
# (fixed) ERROR : Matrix has 33897 rows but found 35891 features (0809 index & t2g, fixed before seurat)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)


Lmatrix <- "../../media/intel6700/Passport/AlicenHonours/work/0909/lambda/seurat_matrix/matrix.mtx.gz"
Lbarcodes <- "../../media/intel6700/Passport/AlicenHonours/work/0909/lambda/seurat_matrix/barcodes.tsv.gz"
Lfeatures <- "../../media/intel6700/Passport/AlicenHonours/work/0909/lambda/seurat_matrix/features.tsv.gz"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)


Umatrix <- "../../media/intel6700/Passport/AlicenHonours/work/0909/untreated/seurat_matrix/matrix.mtx.gz"
Ubarcodes <- "../../media/intel6700/Passport/AlicenHonours/work/0909/untreated/seurat_matrix/barcodes.tsv.gz"
Ufeatures <- "../../media/intel6700/Passport/AlicenHonours/work/0909/untreated/seurat_matrix/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)

