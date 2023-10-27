# Neutrophil  markers 

library(Seurat)
library(ggplot2)
library(reshape2)

# [0] Calculate expression 
# isolate B cells
Neutrophils <- subset(TreatmentAnnotated, celltype == "neutrophils")

markerinfo <- AverageExpression(Neutrophils, 
                                assay = "RNA",
                                features = c(
                                             "OAS1", "OAS2",
                                            "IFIT1", "IFIT2",  
                                            "STAT1", "STAT2",
                                            "MX1",
                                            "THBD","CAMP", "CD93", "VDR"
                                            ),
                                group.by = "treatment",
                                slot = "data",
                                verbose = TRUE)


markerinfo <- as.data.frame(markerinfo)
markerinfo <- t(markerinfo) # get receptors as columns 
rownames(markerinfo) <- sub("^RNA\\.","", rownames(markerinfo))
sample <- rownames(markerinfo)
sample <- as.data.frame(sample)
sample <- data.frame(lapply(sample, function(x) gsub("\\.", " ", x))) # removes dots and replaces with spaces for cell types ONLY

markerinfo <- cbind(sample, markerinfo)

# [1] 

