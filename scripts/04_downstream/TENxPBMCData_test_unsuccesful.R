# Loading the dataset :
library(SingleCellExperiment) # for the dataset
library(TENxPBMCData) # to view the datasets available : args(TENxPBMCData) 
BiocManager::install("SeuratData") #

# Loading the PBMC dataset provided in R by 10X Genomics, using the Read10X function from the Seurat package to load a count matrix from a SingleCellExperiment object.

pbmc3k.data <- TENxPBMCData("pbmc3k") # load dataset into R
counts(pbmc3k.data) # prints the count matrix > empty???

//
# trying to save the count matrix into a file and access it using the tutorial steps :
pbmc3k.count.matrix <- counts(pbmc3k.data)
directory_for_counts <- "~/Alicen/downstream/"
saveRDS(pbmc3k.count.matrix, file.path(directory_for_counts, "pbmc_count_matrix.rds"))

pbmc3k.seurat.object <- CreateSeuratObject(counts = pbmc3k.data) # Create a Seurat object directly from the TENxPBMCData object

# pbmc3k_sce <- SingleCellExperiment(assays = list(counts = pbmc3k), colData = DataFrame(condition = # "pbmc3k")) # create a singlecellexperiment object
# pbmc3k <- Read10X(data.dir = pbmc3k_sce)
