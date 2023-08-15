# This script is to convert ENSEMBL gene names to standard HGNC

# First need to install a mapping between Ensemble gene IDs and HGNC gene names 
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
View(gene_data)
