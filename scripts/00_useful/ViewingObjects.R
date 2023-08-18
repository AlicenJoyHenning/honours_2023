# Doing things with a dataset in R 
# revelation : cache issue : session > clear workspace > close R > open R, reinstall packages

library(readr)
library(dplyr)
library(Matrix)


barcodes <- read_tsv("honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz")
head(barcodes)
length(barcodes)
dim(barcodes)

features <- read_tsv("honours/work/s/features.tsv.gz")
dim(features)
# [1] 69535 2 
# error on reading in the folder containing all these in Read10X() said : 
# Error in dimnamesGets(x, value) : length of Dimnames[[1]] (69536) is not equal to Dim[1] (35639)

# before looking closer at the dimensions, lets rename the columns : 

NewCol <- list('ENSG00000223972.1','DDX11L17') # column heading was an entry, adding the entry to the dataframe
features <- rbind(features, NewCol)
features <- rename(features, EnsembleID = ENSG00000223972.5)
features <- rename(features, HGNC = DDX11L17)

features$EnsembleID <- substr(features$EnsembleID, 1, 15) # removing the '.version' part of the EnsembleID column 

# now we need to address the fact that there are many 'na' value in both columns 

features$HGNC <- ifelse(is.na(features$HGNC), features$EnsembleID, features$HGNC)
features$EnsembleID <- ifelse(is.na(features$EnsembleID), features$HGNC, features$EnsembleID)

# this one wasn't actually necessary 
# features <- features %>%
#   mutate(
#     HGNC = coalesce(HGNC, EnsembleID),
#     EnsembleID = coalesce(EnsembleID, HGNC)
#   )


# Now, looking at the data frame, after row 35638, there are no entries in the ensemble ID column .... just HGNC!?
# check to see if the HGNC names with no ensemble names are duplclates and i can remove them : 

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
# I am going to use this as input for the Read10X function
features <- features[1:35639, ] 

featuresNew <- w

#####

# After addressing problems with the features file, lets look at the problems with the matrix file " 

matrix <- ReadMtx("honours/work/untreated/sm/seurat_matrix/matrix.mtx.gz", "honours/work/untreated/sm/seurat_matrix/barcodes.tsv.gz", "honours/work/untreated/sm/seurat_matrix/features.tsv.gz", skip.feature = 2)
# Error: Matrix has 35639 rows but found 69537 features. Try increasing `skip.feature`. 

matrixx <- writeMM(matrix, "honours/work/untreated/sm/seurat_matrix/newmatrix.mtx.gz")
dim(matrix)
matrix
# 35639 x 5939 sparse Matrix of class "dgCMatrix"


untreated <- CreateSeuratObject( matrix, project='untreated', min.cells=3, min.features=200)
untreated <- CreateSeuratObject(counts=untreated, project='untreated', min.cells=3, min.features=200)


