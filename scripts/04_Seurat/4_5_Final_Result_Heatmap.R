# Load required libraries
library(reshape2)
library(ggplot2)

# Find the average expression for a gene across treatment type:

Bcells <- subset(TreatmentAnnotated, celltype == "B")
AverageExpression(Bcells, 
                  assay = "RNA",
                  features = c("CD79A","CD79B", "CD40", "IFNLR1","IFNAR1", # Receptors
                               
                               
                               
                               ),
                  group.by = "sample",
                  slot = "data",
                  verbose = TRUE)
AverageExpression(TreatmentAnnotated, 
                  assay = "RNA",
                  features = c("CD79A","CD79B" # Receptors
                               
                               
                               
                  ),
                  group.by = "celltype",
                  slot = "data",
                  verbose = TRUE)
