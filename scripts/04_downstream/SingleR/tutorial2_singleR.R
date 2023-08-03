
# 1 install required packages : 
install.packages("SingleR")
install.packages("scran")
install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("SeuratData")
# 2 load required packages : 
library("SingleR")
library("scran")
library("BiocManager")
library(Seurat
BiocManager::install("SingleCellExperiment")
# 3 create and prepare test data : 
test.data <- TENxPBMCData("pbmc4")
# test.counts <- test.data$assays$RNA
seurat_obj <- CreateSeuratObject(counts = test.data$counts, 
                                 assay = "RNA",
                                 meta.data = test.data$meta.data)

# 4 create and prepare reference data : 
HPCAD.ref <- HumanPrimaryCellAtlasData(ensembl = TRUE)
# HPCAD.ref <- CreateSeuratObject(counts = HPCAD.ref$counts, project = "HPCAD.ref")


# 5 perform SingleR annotation : 
test.predictions <- SingleR(test = test.data,
                            assay.type.test = 1,
                            ref = HPCAD.ref,
                            labels = HPCAD.ref$label.main)
# 6 check is predictions contain valid (non-null) annotations 
table(test.predictions$labels)

# 7 create a seurat object for the test data : 



# 7 create UMAP plot with annotations 
DimPlot(test.data, reduction = "umap", group.by = "SingleR.labels")


