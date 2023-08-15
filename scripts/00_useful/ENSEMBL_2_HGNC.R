# This script is to convert ENSEMBL gene names to standard HGNC

##### Install a mapping between Ensemble gene IDs and HGNC gene names ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)

# make a shortcut to access the package : 

edb <- EnsDb.Hsapiens.v86
organism(edb) # making sure it is human 

# create a dataframe from the package with ENSEMBL IDs and HGNC gene names 

gene_data <- select(edb, keys = keys(edb, keytype = "GENEID"), 
                    columns = c("GENEID", "GENENAME"), 
                    keytype = "GENEID", return.only.columns = c("GENEID", "GENENAME"))

# change the names of the columns for clarity : 

colnames(gene_data) <- c("ENSEMBL_Gene_ID", "HGNC_Name")
View(gene_data) # view the wonders 


# save the data frame as a csv 

write.csv(gene_data, file = "honours/gene_data.csv", row.names = FALSE)



##### Create a new feature metadata slot for seurat object (for HGNC) ####

# Convert your mapping data to a named vector : 

ensembl_to_hgnc_map <- setNames(gene_data$HGNC_Name, gene_data$ENSEMBL_Gene_ID)

# Update the feature metadata slot in your Seurat object :

untreated@assays[["RNA"]]@meta.features$gene <- ensembl_to_hgnc_map[rownames(untreated@assays[["RNA"]]@meta.features)]

# Update gene names in the count matrix based on the mapping : 

rownames(untreated@assays[["RNA"]]@counts) <- untreated@assays[["RNA"]]@meta.features$gene

View(untreated) # look 


##### Saving new seurat object #####

saveRDS(untreated, file= "honours/untreated/untreated.rds")
