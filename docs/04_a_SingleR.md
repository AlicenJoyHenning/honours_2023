# Assigning cell types with SingleR 
## Background
SingleR: leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cell independently. Automated annotation method for scRNA-seq data _basically allows you to painlessly assign cell types to scRNA-seq data_



The problem: Without RNA-seq data reference genome > run RNA-seq reads through an assembler to identify sequences > inspect each sequence to determine function (based on sequence motif). This is the case in scRNA-seq > run scRNA-seq cells through clustering protocols > identify cell types.
In scRNA-seq analysis, like a reference genome, the clusters you use to assign cell types act as proxies for the biological states of interest and how accurately you are able to annotate the clusters (assign cell types) is essential. 
Automated cell type annotation methods allow cell type classification by comparing cells ina new dataset to curated reference profiles of known cell types and assigning each new cell to the reference type that its expression profile is most similar to. This allows the scRNA seq analysis to focus on investigating whether or not cell types express changes in adundance or expression across treatments. 

How: 
(1) Calculate Spearman correlation between the expression profile of a cell and the reference sample (?of each experiment, like the untreated, or to the curated dataset cell type?). 
Uses marker genes identified by pairwise comparision between labels in reference data. 
(2) Define per-label score of the correlations across all samples with that label (accounts for differences in sizes of reference samples for each label?)
(3) Repeat score calculation for all labels and find the label with the highest score = SingleR's prediction for the cell 
(4) Improve resolution between closely related labels : scores recalculated with marker genes for subsets of labels, thus focusing on most relevant features, this is iterated until only one label remains. 

## Example
Demonstrating the use of SingleR() on 10X Genomics dataset (Zheng et al. 2017) using Human Primary Cell Atlas dataset (Mabbott et al 2013) as reference. 
This data is loaded using the TENxPBMCData package, specifically designed to provide access to well-known 10X Genomics single-cell RNA-seq datasets of peripheral blood mononuclear cells (PBMCs). 
These datasets are commonly used in single-cell genomics research and analysis. By loading the TENxPBMCData package, you gain access to one or more pre-processed datasets
Similarily, the reference dataset is made available in the package celldex.The r Biocpkg("celldex") package provides convenient access to several cell type reference datasets. Most of these references are derived from bulk RNA-seq or microarray data of cell populations that (hopefully) consist of a pure cell type after sorting and/or culturing. The aim is to provide a common resource for further analysis like cell type annotation of single cell expression data or deconvolution of cell type proportions in bulk expression datasets.

Each dataset contains a log-normalized expression matrix that is intended to be comparable to log-UMI counts from common single-cell protocols [@aran2019reference] or gene length-adjusted values from bulk datasets. By default, gene annotation is returned in terms of gene symbols, but they can be coerced to Ensembl annotation with ensembl=TRUE for more robust cross-referencing across studies.

The datasets available through the TENxPBMCData package are well-known and widely used in the scRNA-seq community. They are preprocessed and formatted in a way that makes them compatible with various downstream analysis tools and workflows.
First install the packages: 
```
# Install the packages :

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TENxPBMCData")
BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("ensembldb") 

# OUTPUT: The downloaded binary packages are in
	C:\Users\alice\AppData\Local\Temp\Rtmp4CNpRY\downloaded_packages

# Load the packages into the R environment :

library(TENxPBMCData) 
library(celldex)
library(SingleR)

# Create objects with the values of the datasets :

new.data <- TENxPBMCData("pbmc4") # using the function TENxPBMCData from the TENxPBMCData package to load a specific dataset
ref.data <- HumanPrimaryCellAtlasData(ensembl=TRUE) # ensures that reference data has Ensembl annotations

# Using SingleR to perform predictions :

predictions <- SingleR(test=new.data, assay.type.test=1, ref=ref.data, labels=ref.data$label.fine, assay.type.test=1)

# create a table to display the results : 

table(predictions$labels)
 


```

```{r}


library(TENxPBMCData)
new.data <- TENxPBMCData("pbmc4k")
