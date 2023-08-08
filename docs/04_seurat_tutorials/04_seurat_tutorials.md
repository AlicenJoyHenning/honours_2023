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
The ```Read10X()``` function from 10X Genomics takes as input the ouptput of 10X cellranger and returns a UMI count matrix: a matrix that represents the number of molecules for each feature (meaning gene > row) that is detected in each cell (column). The count matrix is then used to create a ```Seurat object``` that stores raw and analysis data for a single cell dataset.

```R
# Loading dependencies into your R envirnoment : 
library(dplyr) # part of tidyverse contains functions for working with data frames 
library(Seurat) # package specifically for scRNA-seq
library(patchwork) # customizable layouts for plotting multiple graphs (ggplot2)



## 2) Standard pre-processing workflow 

## 3) Scaling the data 

## 4) Perform linear dimensional reduction 

## 5) Determine the dimensionality of the dataset

## 6) Cluster the cells 

## 7) Run non-linear dimensional reduction (UMAP/tSNE)

## 8) Finding differentially expressed features (cluster biomarkers) 

## 9) Assigning cell type identity to clusters 
