# This code works through singleR & UMAP visualisation of scRNAseq 

# install necessary packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

BiocManager::install("TENxPBMCData")
library(TENxPBMCData)
BiocManager::install("celldex")
library(celldex)
BiocManager::install("SingleR")
library(SingleR)
BiocManager::install("scRNAseq")
library(scRNAseq)
BiocManager::install(c("Biobase", "SummarizedExperiment"))
library(Biobase)
library(SummarizedExperiment)
BiocManager::install("DimPlots")
library(DimPlots)

# create objects associated to the datasets : 

pbmc4_data <- TENxPBMCData("pbmc4")
View(as.data.frame(colData(pbmc4_data)))

HPCAD_ref <- HumanPrimaryCellAtlas(ensembl=TRUE)
View(as.data.frame(colData(HPCAD_ref)))

# execute SingleR to get cell type predictions 

predictions <- SingleR(test = pbmc4_data, ref = HPCAD_ref, labels = HPCAD.ref$labels.main, assay.type.test = 1)
View(as.data.frame(colData(predictions)))

# visualize the cell predictions using the labels assigned by SingleR 

pbmc4_data$SingleR.labels <- predictions$labels[match(rownames(pbmc4_data), rownames(predictions))]
DimPlot(pbmc4_data, reduction = "umap", group.by = "SingleR.labels")
