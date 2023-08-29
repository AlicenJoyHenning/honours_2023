# Differential Gene Expression Analysis for Single Cell RNA sequencing 
## BACKGROUND 
**Differential Gene Expression**: <br> <br>
RNA sequencing provides huge amounts of information on gene expression, often in the context of various conditions. This leads to the concept of _differential_ gene expression (*D*GE) meaning the expression of a gene is different when compared across two or more conditions. The same is true for single cell RNA sequencing where DGE is a process that helps us find out which genes are working differently depending on the cell-type specific conditions we're looking at. In both bulk and single cell RNA sequencing, differentially expressed genes are identified according to statistical signficance of the expressed genes. This information is then best summarised by generating figures. However, the tools used to do so for bulk and single cell RNA sequencing differ.  

## Bulk RNA sequencing DGEA

In bulk RNA-seq, *Cuffdiff, edgeR, and DESeq2* are the most widely used tools to generate figures and have been combined into one ppackage : Visualisation of Differential Gene Expression Results using R (ViDGER)[1]. The functions included in this package are described below and fall into one of three tiers: 
- Tier I : visualise raw/normalised expression levels and display trends 
- Tier 2 :  utilise statistical tests to provide information at the specfic gene level. 

| Function | Tier | Role | Plots |
|----------|------|------|-------|
| Treatment Distributions| I | Detect bias introduced by treatment: if overall distributions are not similar | boxplot / violin |
| Comparision of expression levels | I | Compare gene expression across 2 conditions where points (genes) that occur above/below the diagonal indicate that the gene has higher/lower expression for the y axis condition relative to the x axis condition | *scatter plot (cond 1 vs cond 2) | 
| Number of DGEs | I | Shows the number of pairwise DEGs for all treatment conditions indicating which comparisions are more dissimilar | histogram/heatmap |
| Fold-change vs normlaised mean counts | II | FC vs NMC between two conditions with log FC on y axis and MC on x axis | scatter plot |
| P-value vs fold-change | II | Statistical significance of the difference in gene expression between 2 genes relative to the magnitufe od the difference in gene expression for ever gene in the comparision: higher up on y axis meaning smaller p value, overall wider dispersion indicates two treatment groups have higher levels of difference regarding gene expression | Volcano Plot |
| Four-way plot | II  | relative comparison of fold change: two distinct treatment groups compared relative to a control | segmented scatter plot | 

*must include the normalized expression values (FPKM/CPM/base-10 log)

## Single Cell RNA sequencing DGEA

In single cell cases, DGEA aims to identify genes that are differentially expressed between cell subpopulations and between experimental conditions. Some methods to find DEG include DESingle [2] and seurat functions (FindMArkers()). 

Seurat's FindMarkers function is part of the Seurat package, which is widely used for scRNA-seq analysis. This function uses a negative binomial (NB) regression model to identify differentially expressed genes between two groups (condition and control). It estimates dispersion and performs statistical tests to assess significance. It's worth noting that Seurat is a versatile package with a wide range of analysis capabilities, not solely focused on differential expression.

DEsingleR, on the other hand, is specifically designed for single-sample differential expression analysis, which makes it relevant for scenarios where replicates are lacking. DEsingleR estimates biological and technical variability from single-cell data and applies a likelihood ratio test to identify differentially expressed genes between two conditions. It's designed to handle the unique challenges of single-sample comparisons in scRNA-seq data.

However, because the output of the seurat differential gene expression functions are compatiable with the downstream gene set enrichment analysis, it was chosen for this study. 





## INTERPRETATION 



## REFERENCES 
[1] McDermaid, A., Monier, B., Zhao, J., Liu, B. and Ma, Q., 2019. Interpretation of differential gene expression results of RNA-seq data: review and integration. Briefings in bioinformatics, 20(6), pp.2044-2054.
[2] Ma, Y., Sun, S., Shang, X., Keller, E.T., Chen, M. and Zhou, X., 2020. Integrative differential expression and gene set enrichment analysis using summary statistics for scRNA-seq studies. Nature communications, 11(1), p.1585.
