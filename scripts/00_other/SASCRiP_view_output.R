# VIEWING SASCRiP OUTPUT FILES

# revelation : cache issue : session > clear workspace > close R > open R, reinstall packages

##### [1] Viewing features ####
library(readr)
library(dplyr)
library(Matrix)
library(stringr)
library(Seurat)


barcodes <- read_tsv("honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz")
head(barcodes)
length(barcodes)
dim(barcodes)

UntreatedFeatures <- read_tsv("honours/work/untreated_my_index/AdjustedFeatures.tsv.gz", col_names = FALSE)
dim(UntreatedFeatures)
# [1] 69535 2 
# error on reading in the folder containing all these in Read10X() said : 
# Error in dimnamesGets(x, value) : length of Dimnames[[1]] (69536) is not equal to Dim[1] (35639)

# before looking closer at the dimensions, lets rename the columns : 

# NewCol <- list('ENSG00000223972.1','DDX11L17') # column heading was an entry, adding the entry to the dataframe
# UntreatedFeatures <- rbind(NewCol, UntreatedFeatures)
UntreatedFeatures <- rename(UntreatedFeatures, Ensembl_ID = X1)
UntreatedFeatures <- rename(UntreatedFeatures, HGNC = X2)
colnames(UntreatedFeatures)


UntreatedFeatures$Ensembl_ID <- substr(UntreatedFeatures$X1, 1, 15) # removing the '.version' part of the EnsembleID column 
UntreatedFeatures$HGNC <- UntreatedFeatures$X2
UntreatedFeatures <- UntreatedFeatures %>% select(-1)
UntreatedFeatures <- UntreatedFeatures %>% select(-1)
# now we need to address the fact that there are many 'na' value in both columns 

UntreatedFeatures$HGNC <- ifelse(is.na(UntreatedFeatures$HGNC), UntreatedFeatures$Ensembl_ID, UntreatedFeatures$HGNC)
# features$EnsembleID <- ifelse(is.na(features$EnsembleID), features$HGNC, features$EnsembleID)

# this one wasn't actually necessary 
# features <- features %>%
#   mutate(
#     HGNC = coalesce(HGNC, EnsembleID),
#     EnsembleID = coalesce(EnsembleID, HGNC)
#   )

# Now, looking at the data frame, after row 35638, there are no entries in the ensemble ID column .... just HGNC!?
# check to see if the HGNC names with no ensemble names are duplclates and i can remove them : 

# save new features file 
UntreatedFeatures <- write_tsv(UntreatedFeatures, file = "honours/work/untreated_my_index/NewFeatures.tsv.gz")
UntreatedFeatures <- write_tsv(UntreatedFeatures, file = "honours/work/untreated/seurat_matrix/UnduplicatedFeatures.tsv.gz")

##### [2] Duplication Problem ####

head(duplicated(features$HGNC)) 
# there is a true!! sad. 
head(which(duplicated(features$HGNC)))
# this finds the position of indices of all duplicated elements 
duplicated <- which(duplicated(features$HGNC))
length(duplicated)      
# 31990 !!! hectic
# 26112 after adjustments of replacing NAs ... duh cmon Alicen
duplicated <- which(duplicated(features$EnsembleID))
length(duplicated)
# 17 .... what? surely it should be the same as the repeated HGNC... or not i guess if hgnc was used to replace an ensemble ID 

 
# Actually not going to remove the duplicates becasue i don't know if it searches for duplicates across columns, also i don't know if it keeps one of the values : which i need it to! So atm i'm okay with keeping the duplicates.
# featuresNoReplicates <- features[!duplicated(features$HGNC), ] # subset of the features object that excludes duplicated rows based on the HGNC column
# # 37546 BUT still has empty ensemble ID column for many genes 

#featuresNoReplicates <- features %>% distinct(HGNC)
#dim(featuresNoReplicates) # ouput is a subset of distinct entries OF THAT column, not the entire dataset 

# I am going to use this as input for the Read10X function to see if it works : 
features <- features[1:35639, ] 
# it doesn't work into the Read10X function, but it works when i create a matrix and and use that directly to make a seurat object 
# but clusters are weird so 

##### [3] Subset based on features from other datasets ####

UntreatedFeatures <- read_tsv("honours/work/untreated/seurat_matrix_unsuccessful/features.tsv.gz")
dim(UntreatedFeatures)
# [1] 69535 2 
AlphaFeatures <- read_tsv("honours/work/ifnalpha/seurat_matrix/features.tsv.gz")
dim(AlphaFeatures)
# Rows: 33896 Columns: 2   
LambdaFeatures <- read_tsv("honours/work/ifnlambda/seurat_matrix/features.tsv.gz")
# Rows: 33896 Columns: 2   

# clearly (due to the dimensions being almost double) there is something wrong with the features of untreated. 
# plan : subset the UnTreated Features using the features from alpha and lambda 

sum(str_detect(UntreatedFeatures$HGNC,'PRDM16'))
# [1] 2 > - we this gene is repeated in untreated features. 
sum(str_detect(UntreatedFeatures$HGNC,'RPS27P9'))
# [1] 1
# This tells me that removing duplicates might be a good idea : 

duplicated <- which(duplicated(UntreatedFeatures$HGNC))
length(duplicated) 
# [1] 26112

UntreatedFeatures <- UntreatedFeatures[!duplicated(UntreatedFeatures$HGNC), ] 
dim(UntreatedFeatures)
# 43424     2  better still significantly more than alpha and lambda features 


sum(str_detect(UntreatedFeatures$HGNC,'PRDM16'))
# [1] 1 Good 


UntreatedFeatures <- rename(UntreatedFeatures, Ensembl_ID = EnsembleID)

UntreatedFeatures <- write_tsv(UntreatedFeatures, file = "honours/work/untreated/seurat_matrix/UnduplicatedFeatures.tsv.gz")
Test <- merge(UntreatedFeatures, AlphaFeatures, by = "Ensembl_ID", all.x = TRUE)
# output shows that none of the HGNC names for the same Ensemble id are the same 

Test <- Test %>% select(-HGNC.x)
Test <- rename(Test, HGNC = HGNC.y)
Test <- write_tsv(Test, file = "honours/work/untreated/seurat_matrix_unsuccessful/HGNCfeatures.tsv.gz")

UntreatedTest <- read_tsv("honours/work/untreated/sm_new_index_HGNC/features.tsv.gz")
dim(UntreatedFeatures)

##### [4] Making Matrix #####

# After addressing problems with the features file, lets look at the problems with the matrix file " 

matrix <- ReadMtx("honours/work/untreated/seurat_matrix_unsuccessful/matrix.mtx.gz", "honours/work/untreated/seurat_matrix_unsuccessful/barcodes.tsv.gz", "honours/work/untreated/seurat_matrix_unsuccessful/features.tsv.gz")
# Error: Matrix has 35639 rows but found 69537 features. Try increasing `skip.feature`. 
# Error : Error: Matrix has 35639 rows but found 43425 features, 7786 features need to be removed 

# (Edit: can't) extract the features from matrix file and use them to subset the features file 

UntreatedMatrix <- readMM("honours/work/untreated/sm/seurat_matrix/matrix.mtx.gz")
dim(UntreatedMatrix)
UntreatedGenes <- rownames(UntreatedMatrix) # empty 

AlphaMatrix <- readMM("honours/work/ifnalpha/seurat_matrix/matrix.mtx.gz")
AlphaGenes <- rownames(AlphaMatrix) # also empty 




# matrixx <- writeMM(matrix, "honours/work/untreated/sm/seurat_matrix/newmatrix.mtx.gz")
# dim(matrix)
# matrix
# # 35639 x 5939 sparse Matrix of class "dgCMatrix"


untreated <- CreateSeuratObject(matrix, project='untreated', min.cells=3, min.features=200)
untreated <- CreateSeuratObject(counts=untreated, project='untreated', min.cells=3, min.features=200)

UntreatedMatrix <- ReadMtx("honours/work/untreated/sm/seurat_matrix/matrix.mtx.gz", "honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz", "honours/work/ifnalpha/seurat_matrix/features.tsv.gz", skip.feature = 33900)

##### [5] Again #####

UntreatedFeatures <- read_tsv("honours/work/s/features.tsv.gz")
dim(UntreatedFeatures)
# 35640 

UntreatedMatrix <- readMM("honours/work/s/matrix.mtx.gz")
dim(UntreatedMatrix)
# 35 639 5939 

matrix <- ReadMtx("honours/work/s/matrix.mtx.gz", "honours/work/s/barcodes.tsv.gz", "honours/work/s/features.tsv.gz", skip.feature = 2)
untreated <- CreateSeuratObject(matrix, project='untreated', min.cells=3, min.features=200)


