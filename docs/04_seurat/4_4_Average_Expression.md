## AverageExpression() Seurat function

Using this function, you can find the average expression according to any existing or newly created meta data column for any collection or subset of cells. Here, I am looking at a B cell cluster (subsetted from the original Seurat object) and finding the expression according to treatment condition (stored in the metadata column) : 

``` 
# Subsetting Seurat Object to isolate B cell : 
Bcells <- subset(TreatmentAnnotated, celltype == "B")

# Finding the avg. exp in B cells according to treatment : 
AverageExpression(Bcells, 
                  assay = "RNA",
                  features = c("CD79A","CD79B"),
                  group.by = "sample",
                  slot = "data",
                  verbose = TRUE)
```
OUTPUT: 
```
            alpha     lambda   untreated
CD79A  8.02030489 13.3625480 13.30160284
CD79B  0.96551132  4.0050962  3.18476000
```

Another example would be to look at the expression of the same genes but in accordance to cell types from the original Seurat object : 

```
AverageExpression(TreatmentAnnotated, 
                  assay = "RNA",
                  features = c("CD79A","CD79B" ),
                  group.by = "celltype",
                  slot = "data",
                  verbose = TRUE)
```
OUTPUT:  
```
       monocytes naive CD4 T neutrophils   T helper naive CD8 T        B cytotoxic T        mDCs
CD79A 0.03379605  0.01850732  0.04400375 0.02404872  0.06299077 9.514192  0.03025184 0.033897478
CD79B 0.02483115  0.02940449  0.00852869 0.04934424  0.02094900 1.682108  0.04425225 0.009832255

              NKT      Tregs         NK platelets unknown       DCs        Tcm       pDCs
CD79A 0.004120578 0.02564899 0.01899818 0.0237351       0 0.3136440 0.00000000 0.03510176
CD79B 0.059715230 0.04655470 0.15528257 0.0000000       0 0.2423257 0.01722706 0.00000000
```
