
# Assigning cell type identity to clusters 

## (1) Canonical markers 
We can use canonical markers to match the unbiased clustering to known cell types. This requires the differentially expressed features of each cluster to be founda and compared to these.
I found this method difficult. 

| Markers | Cell Type |
|----------|----------|
| IL7R, CCR7 | Naive CD4+ T   |
| CD14, LYZ	CD14+   | Monocytes   |
| IL7R, S100A4	| Memory CD4+  |
|	MS4A1	| B cells |
|	CD8A	| CD8+ T |
|	FCGR3A, MS4A7	FCGR3A+ | Monocytes |
|	GNLY, NKG7	| NK | 
|	FCER1A, CST3 |	DC |
|	PPBP	| Platelet | 

```R
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet", "Other1", "Other2")
names(new.cluster.ids) <- levels(alpha)
alpha <- RenameIdents(alpha, new.cluster.ids)
DimPlot(alpha, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```
![image](https://github.com/AlicenJoyHenning/honours_2023/blob/main/plots/alpha_clusters_1.jpg)



## (2) scAnnotatR 
As an alternative, the scAnnotatR package performs automatic cluster annotation based on pretrained models. The problem I am having with this at the moment is that it is seems to assign cell types based on the presence of all marker genes, meaning that if a cluster doesn not have all the marker genes for a cell type it will not be classified as such. However, this resulted in no clusters being annotated because none of them displayed all marker genes. 

```R
BiocManager::install("devtools")
BiocManager::install("scAnnotatR")

library(scAnnotatR)

default_models <- load_models("default") # loading the pre-trained models 
names(default_models)
# [1] "B cells"           "Plasma cells"      "NK"                "CD16 NK"           "CD56 NK"           "T cells"          
# [7] "CD4 T cells"       "CD8 T cells"       "Treg"              "NKT"               "ILC"               "Monocytes"        
# [13] "CD14 Mono"         "CD16 Mono"         "DC"                "pDC"               "Endothelial cells" "LEC"              
# [19] "VEC"               "Platelets"         "RBC"               "Melanocyte"        "Schwann cells"     "Pericytes"        
# [25] "Mast cells"        "Keratinocytes"     "alpha"             "beta"              "delta"             "gamma"            
# [31] "acinar"            "ductal"            "Fibroblasts"

# takes as input a seurat object : 
is(alpha, "Seurat")
# [1] TRUE

# To launch cell type identification, we simply call the `classify_cells`function : 
alpha.scannotatR <- classify_cells(classify_obj = alpha, 
                             assay = 'RNA', slot = 'counts',
                             cell_types = 'all', 
                             path_to_models = 'default')

# ERROR : All genes from "" classifier model must be present in the dataset to perform classification. Classification of B "" skipped.

```


## (3) SingleR 



