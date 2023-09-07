# Viewing the SASCRiP output 

library(dplyr)

##### Read in files #####
UntreatedFeatures <- read_tsv("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz", col_names = FALSE)
dim(UntreatedFeatures)
# Rows: 35639 Columns: 2  

AlphaFeatures <- read_tsv("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/features.tsv.gz", col_names = FALSE)
dim(AlphaFeatures)
# [1] 35639     2

LambdaFeatures <- read_tsv("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/features.tsv.gz", col_names = FALSE)
dim(LambdaFeatures)
# Rows: 35639 Columns: 2 
##### Fix column names #####

UntreatedFeatures$Ensembl_ID <- substr(UntreatedFeatures$X1, 1, 15)
UntreatedFeatures$HGNC <- UntreatedFeatures$X2
UntreatedFeatures <- UntreatedFeatures %>% select(-1) # times 2 
UntreatedFeatures <- write_tsv(UntreatedFeatures, file = "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz")

AlphaFeatures$Ensembl_ID <- substr(AlphaFeatures$X1, 1, 15)
AlphaFeatures$HGNC <- AlphaFeatures$X2
AlphaFeatures <- AlphaFeatures %>% select(-1) # times 2 
AlphaFeatures <- write_tsv(AlphaFeatures, file = "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz")

LambdaFeatures$Ensembl_ID <- substr(LambdaFeatures$X1, 1, 15)
LambdaFeatures$HGNC <- LambdaFeatures$X2
LambdaFeatures <- LambdaFeatures %>% select(-1) # times 2 
LambdaFeatures <- write_tsv(LambdaFeatures, file = "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz")

