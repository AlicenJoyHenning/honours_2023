# Alpha Dataset Processing in R with Seurat 

## Preprocessing Stages  
After going through the first two functions of SASCRiP, the output folder containing barcodes, features, and matrix files were used to go through a manual preprocessing workflow offered by Seurat.

```R

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(grid)

getwd()

alpha.data <-  Read10X(data.dir = "honours/ifnalpha/seurat_matrix/")
alpha <- CreateSeuratObject(counts=alpha.data, project='ifnalpha', min.cells=3, min.features=200)
# output : Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')

alpha
# An object of class Seurat 
# 17692 features across 6158 samples within 1 assay 
# Active assay: RNA (17692 features, 0 variable features)

alpha.data[1:3, 1:3]
# 3 x 3 sparse Matrix of class "dgCMatrix"
# AAACCCAAGACTCATC AAACCCACAACGTATC AAACCCACAATTGCTG
# HGNC                  .                .                .
# PRDM16                .                .                .
# PEX10       
```

Next, the dataset must be altered to remove low quality cells that are determined by the number of features expressed and the mitochondrial percentage

```R

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")

alpha <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
```
After quality control, the dataset is normalized : 

```R

# 1 : actual normalization 
alpha <- NormalizeData(alpha, normalization.method = "LogNormalize", scale.factor = 10000)
# Performing log-normalization
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|

# 2 : feature selection 
alpha <- FindVariableFeatures(alpha, selection.method = "vst", nfeatures = 2000)
# Calculating gene variances
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Calculating feature variances of standardized and clipped values
#    0%   10   20   30   40   50   60   70   80   90   100%
#      [----|----|----|----|----|----|----|----|----|----|
#         **************************************************|

# 3 : scaling the data 
all.genes <- rownames(alpha)
alpha <- ScaleData(alpha, features = all.genes)
# Centering and scaling data matrix
# |==========================================================================| 100%
# The results of this are stored in pbmc[["RNA"]]@scale.data

```

## Linear dimensional reduction 

Next perform PCA on the scaled data. By default, only the previously determined variable features are used as input

```R alpha <- RunPCA(alpha, features = VariableFeatures(object = alpha))```

This generates the output: 
```R
PC_ 1 
Positive:  RPL23A, RPL3, RPL13A, RPL26, RPS3A, RPL18A, RPS23, RPL7A, RPS18, RPL35 
ENSG00000237550, RPSA, RPL10A, RPL5, RPL10, RPS4X, EEF1B2, RPL6, RPS7, RPS8 
RPS5, RPL18, RPS20, RPL37A, RPS15A, RPL12, RPL24, RPL35A, RPS6, RPS3 
Negative:  IL1RN, S100A8, FTH1, FFAR2, MNDA, ACSL1, APOBEC3A, TYROBP, NCF1C, NINJ1 
GLUL, SNX10, IFITM3, FCER1G, SPI1, NCF1, TYMP, C15orf48, IER3, TLR2 
CCL4, CLEC7A, CCL4L2, LST1, S100A11, NCF2, VNN2, SAMSN1, ISG15, BASP1 
PC_ 2 
Positive:  LEF1, LDHB, IL7R, BCL11B, IKZF1, TCF7, GPR155, ADTRP, ITGA6, CAMK4 
GIMAP4, CMTM8, LEPROTL1, C12orf57, MAL, PRMT2, PDCD4, CD3G, CCR7, TNFRSF25 
MYC, SCML1, LTB, STING1, GIMAP5, FHIT, FAM184A, GIMAP6, PASK, ABLIM1 
Negative:  HLA-DPA1, HLA-DRA, CD74, MYOF, KYNU, HLA-DQB1, LILRB1, CD86, RGL1, HLA-DQA1 
LGALS1, ADAM28, EPB41L3, CPVL, CYP1B1, OLR1, ANXA2, LACC1, EMP1, LY86 
BANK1, MS4A1, TCF4, MS4A7, RIN2, CCL2, CCL7, SWAP70, ENG, FGD2 
PC_ 3 
Positive:  MS4A1, IGHM, BANK1, BCL11A, IGHD, IGKC, HLA-DQA1, HLA-DPA1, ADAM28, HLA-DRA 
TCF4, IGKJ5, CCSER1, CD19, PLEKHG7, RAB30, BCL2, MEF2C, NAPSB, EAF2 
HLA-DQB1, CCR7, IGLC2, CD79B, HLA-DPB1, STAP1, CD24, ACSM3, CD74, C7orf50 
Negative:  NKG7, PRF1, GNLY, ANXA1, GZMB, CLIC3, SAMD3, MYBL1, METRNL, SYTL2 
CD8A, FYN, SH2D2A, C1orf21, KLRD1, IL18RAP, GZMM, CD96, PYHIN1, PLAAT4 
SYTL3, AUTS2, EOMES, KLRF1, IQGAP2, PIK3R1, LAG3, MAF, TNFRSF18, TRGC2 
PC_ 4 
Positive:  NKG7, PRF1, GNLY, GZMB, CLIC3, SAMD3, KLRD1, MS4A1, KLRF1, BANK1 
IGHM, BCL11A, IL18RAP, MYBL1, CD38, C1orf21, IGHD, SH2D2A, EOMES, ADAM28 
CD8A, IGKC, STAP1, TNFRSF18, SYTL2, CCSER1, PYHIN1, FGFBP2, GZMM, AUTS2 
Negative:  MYOF, EPB41L3, RGL1, CYP1B1, CPVL, EMP1, OLR1, CCL2, CCL7, LACC1 
CD86, CXCL11, ENG, MS4A7, PLA2G7, CTSB, P2RY6, LYZ, MSR1, LILRB4 
TGFBI, LILRB1, LEF1, SRC, LGALS1, CST3, KYNU, RIN2, THBS1, RGS10 
PC_ 5 
Positive:  ENSG00000284874, NRGN, ACRBP, TSC22D1, DAB2, PTCRA, TUBA4A, CMTM5, MMD, TMEM40 
C2orf88, ENKUR, CLU, PGRMC1, CTTN, H2BC11, MPIG6B, MTURN, SPARC, PDE5A 
GTF3C2, MAP3K7CL, H2AJ, TNFSF4, CLDN5, BEX3, H2AC6, CTSA, SSX2IP, TSPAN33 
Negative:  ISG15, MYOF, RGL1, LGALS1, S100A11, APOBEC3A, LILRB1, TYROBP, OLR1, MS4A7 
LACC1, RIN2, LILRB4, CPVL, KYNU, EMP1, IL1RN, ANXA5, CTSB, S100A4 
IFITM3, ENG, CYP1B1, MT2A, CCL4, PLA2G7, S100A6, GCH1, IL4I1, CD86 
```
This code performs Principal Component Analysis (PCA) on a Seurat object named alpha.
alpha: This is the Seurat object on which you want to perform PCA.
RunPCA(): This is a function provided by the Seurat package that performs PCA. PCA is a dimensionality reduction technique that is commonly used to reduce the complexity of high-dimensional data while retaining as much information as possible.
features = VariableFeatures(object = alpha): This part specifies the features (genes) to be used for PCA. It uses the VariableFeatures() function to select the highly variable features from the Seurat object alpha. Highly variable features are those that exhibit the most biological variation across cells, making them important for capturing cell-to-cell differences.

Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. 

![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/9f67f1ab-1766-48cf-b795-b68ee890a7b6)

## Determine dimensionality of the dataset 

In scRNA-seq data, there can be a lot of noise (variability not related to biological differences) in the gene expression measurements. To make sense of the data and find meaningful patterns, it's important to reduce this noise. One way to do this is by using Principal Component Analysis (PCA), which finds the most important combinations of genes that capture the major sources of variation in the data.

```R
alpha <- JackStraw(alpha, num.replicate = 100)
alpha <- ScoreJackStraw(alpha, dims = 1:20)
```
The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features with low p-values (solid curve above the dashed line). JackStraw Procedure: The JackStraw procedure involves randomly shuffling the data a bit (permuting) and then re-running PCA. By doing this multiple times, you create a "null distribution" of how the scores might look if there's no real signal. PCs that significantly deviate from this null distribution are considered significant.

It shows the distribution of p-values for each PC. A p-value is a measure of how likely it is to observe a certain value by random chance. If a PC's p-value is very low (solid curve above the dashed line), it means that the observed score is unlikely to have occurred by random chance alone. These significant PCs are the ones you should pay attention to.

By identifying the significant PCs, you're ensuring that the features (genes) you're using for clustering cells are meaningful and not dominated by noise. This step helps you focus on the most important sources of variation and makes your subsequent analyses, like clustering, more accurate and biologically relevant.

![image](http://127.0.0.1:18611/graphics/6da99052-2642-4a4b-b6e2-a5ecc551e36c.png)






