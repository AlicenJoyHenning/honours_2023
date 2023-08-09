# Seurat -Guided Clustering Tutorial 
# Contents
- <span style="color: red;">[(1) Setup the Seurat Object](#section-1)</span>
- <span style="color: orange;">[(2) Standard pre-processing workflow](#section-2)</span>
- <span style="color: yellow;">[(3) Scaling the data](#section-3)</span>
- <span style="color: green;">[(4) Perform linear dimensional reduction](#section-4)</span>
- <span style="color: blue;">[(5) Determine the dimensionality of the dataset](#section-5)</span>
- <span style="color: indigo;">[(6) Cluster the cells](#section-6)</span>
- <span style="color: violet;">[(7) Run non-linear dimensional reduction (UMAP/tSNE)](#section-7)</span>
- <span style="color: purple;">[(8) Finding differentially expressed features (cluster biomarkers)](#section-8)</span>
- <span style="color: pink;">[(9) Assigning cell type identity to clusters](#section-9)</span>

## 1) Setup the Seurat Object 
The raw data was downloaded and extracted easily on a Linux machine but in Windows had to be extracted using the command ```tar -xvzf C:\PATH\TO\FILE\FILE-NAME.tar.gz -C C:\PATH\TO\FOLDER\EXTRACTION```.

The ```Read10X()``` function from 10X Genomics takes as input the ouptput of 10X cellranger and returns a UMI count matrix: a matrix that represents the number of molecules for each feature (meaning gene > row) that is detected in each cell (column). The count matrix is then used to create a ```Seurat object``` that stores raw and analysis data for a single cell dataset.

```R

# Loading dependencies into your R envirnoment : 

library(dplyr) # part of tidyverse contains functions for working with data frames 
library(Seurat) # package specifically for scRNA-seq
library(patchwork) # customizable layouts for plotting multiple graphs (ggplot2)

#  Load the dataset from saved location :

pbmc.data <-  Read10X(data.dir = "honours_2023/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(
counts=pbmc.data, # specifies the matrix of raw counts representing the gene expression data for each cell (row > gene, column > cell).
project='pbmc3k', # assign a project name or identifier to the Seurat objecT
min.cells=3, # specifies the minimum number of cells a gene must be detected in to be considered for inclusion in the analysis (FILTERS OUT GENES)
min.features=200 # sets the minimum number of features (genes) a cell must express to be considered for inclusion in the analysis (FILTERS OUT CELLS)
)

pbmc
# An object of class Seurat 
# 13714 features across 2700 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# viewing the data as a sparse matrix (dots save space)
# CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
# TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
# MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .

pbmc.data[1:5, 1]
# viewing expression of first 5 genes in the first cell (when only one cell, no longer sparse matrix)
# MIR1302-10      FAM138A        OR4F5 RP11-34P13.7 RP11-34P13.8 
# 0            0            0            0            0

pbmc.data[c("CD3D"), 1:5]
# viewing the expression of one gene in first 5 cells : 
# AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 AAACCGTGTATGCG-1 
# 4                0               10                0                0

```

## 2) Standard pre-processing workflow 

To explore the data more and answer questions like which genes are most highly expressed across all cell types, we need to first preprocess the data. 

```R

```


## 3) Scaling the data 

## 4) Perform linear dimensional reduction 

## 5) Determine the dimensionality of the dataset

## 6) Cluster the cells 

## 7) Run non-linear dimensional reduction (UMAP/tSNE)

## 8) Finding differentially expressed features (cluster biomarkers) 

## 9) Assigning cell type identity to clusters 
