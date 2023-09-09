# EDIT T2G FILE FOR SASCRiP 

##### [1] Calling in the data #####
# Load dependencies :
install.packages("readr")
library(readr)
library(Seurat)

# Data stored in the tsv file using read_tsv package and assigning the data to a 'variable'(check name)

features <- read_tsv("D:/AlicenHonours/kallisto_index/0809version/t2g.txt", col_names = FALSE)
darisia_features <- read_tsv("D:/AlicenHonours/darisia_index/transcripts_to_genes.txt", col_names = FALSE)

##### [2] Edit column names #####

# check names of columns : 
colnames(features)

# change the names of the headers : 
colnames(features) <- c("Transcript_ID","Gene_ID", "HGNC")

# check if it worked :  
colnames(features)

# UPDATE : we don't remove version names because they are helpful for assigning the correct gene names (changes do happen & this allows you to keep track)
# Remove the version number from the ensembl ID gene names before using them for the HGNC names 
# features$Ensembl_ID <- substring(features$Ensembl_ID, 1, 15)

# Replace empty 'na' rows in the HGNC column with gene ID 
# if else statement to change those 'NA' values in the GGNC column with corresponding ENSG name 
# ifelse(test (what condition),yes (if condition holds do this),no (if condition doesn't hold do this))
features$HGNC <- ifelse(is.na(features$HGNC),features$Gene_ID, features$HGNC)

##### [3] save the updated version of the tsv file #####

write_tsv(features, "D:/AlicenHonours/kallisto_index/0909version/corrected_t2g.tsv.gz")

# check updated_features file : 
updated_features <- read_tsv("updated_features.tsv.gz")

##### [4] To do the opposite #####
# From the column where NA was replaced with ENSG names, I want the NA to be returned : 
  
  features$HGNC <- ifelse(grepl("^ENSG", features$HGNC), NA, features$HGNC)