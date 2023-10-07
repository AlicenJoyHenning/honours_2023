# DIFFERENTIAL EXPRESSION ANALYSIS (Seurat) 

##### [1.1] Dependencies & load object #####
BiocManager::install("clusterProfiler", version = "3.14")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("gage")
BiocManager::install("gageData")

library(clusterProfiler)
library(pathview)
library(DOSE)
library(biomaRt) # for entrezgenes (not 'mart')
library(Seurat)
library(SeuratObject)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(patchwork)
library(openxlsx)


TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
Idents(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated) <- "RNA"

##### [1.2] View the annotated clusters #####
# View the current annotations for confirmation : 
Idents(TreatmentAnnotated) 
# Defined colour palette adjusted for cluster annotation : 
# Visualise UMAP plot with annotations : 
all <- DimPlot(TreatmentAnnotated, reduction = "umap", pt.size = 1.5, label = TRUE, label.color = "white", label.size = 6, label.box = TRUE, repel = TRUE, cols = palette.b)

# Quick access feature plot to check for gene presence : 
features <- c("CXCR1", "CXCR2", "CD4", "KLF2", "IL7R", "CD8B", "CD8A", "BCL11A", "CD40", "NKG7", "CTSW", "FOXP3", "RTKN2")
f <- DotPlot(object = treatment, 
        features = features,
        cols = c("grey", "#6ab5ba")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

# Patchwork plot with total UMAP with feature plots : 
all | f

##### [2.1] DE with FindMarkers() on fine cell types (15) ######
# Ensure the default assay is RNA not integrated : 
DefaultAssay(TreatmentAnnotated)
DefaultAssay(TreatmentAnnotated) <- "RNA"

# Creating metadata column to hold cell type AND stimulation information (as opposed to just cell type) : 
TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$treatment, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the idents to that column 

# Now we can easily refer to a cell & treatment type, to see the options : 
levels(TreatmentAnnotated) # list the options to perform DE on :

# Use FindMarkers() to find the genes that are different between stimulated and untreated cell types

# MYELOID CELL TYPES : 
# 1. Monocytes 
M1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "monocytes_alpha", ident.2 = "monocytes_untreated", sep = "_")
M1AlphaResponse$gene <- rownames(M1AlphaResponse) # puts gene names into a column
M1AlphaResponse <- filter(M1AlphaResponse, p_val_adj < 0.05)
down <- filter(M1AlphaResponse, avg_log2FC < 0)
up <- filter(M1AlphaResponse, avg_log2FC >= 0)
write.table(M1AlphaResponse$gene, file = file_path, col.names = FALSE, row.names = FALSE)
# saveRDS(M1AlphaResponse, "honours/results/FinalIndex/DEAnalysis/monocytesAlphaResponse.rds")

M1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "monocytes_lambda", ident.2 = "monocytes_untreated", sep = "_") 
M1LambdaResponse$gene <- rownames(M1LambdaResponse) 
M1LambdaResponse <- filter(M1LambdaResponse, p_val_adj < 0.05)
down <- filter(M1LambdaResponse, avg_log2FC < 0)
up <- filter(M1LambdaResponse, avg_log2FC >= 0)
# saveRDS(M1LambdaResponse, "honours/results/FinalIndex/DEAnalysis/monocytesLambdaResponse.rds")


# 2. Neutrophils 
M2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_alpha", ident.2 = "neutrophils_untreated", sep = "_") 
M2AlphaResponse$gene <- rownames(M2AlphaResponse)
M2AlphaResponse <- filter(M2AlphaResponse, p_val_adj < 0.05)
down <- filter(M2AlphaResponse, avg_log2FC < 0)
up <- filter(M2AlphaResponse, avg_log2FC >= 0)
saveRDS(M2AlphaResponse, "honours/results/FinalIndex/DEAnalysis/NeutrophilsAlphaResponse.rds")

M2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_lambda", ident.2 = "neutrophils_untreated", sep = "_") 
M2LambdaResponse$gene <- rownames(M2LambdaResponse) 
M2LambdaResponse <- filter(M2LambdaResponse, p_val_adj < 0.05)
down <- filter(M2LambdaResponse, avg_log2FC < 0)
up <- filter(M2LambdaResponse, avg_log2FC >= 0)
saveRDS(M2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/NeutrophilsLambdaResponse.rds")

# DENDRITIC CELL TYPES 
# 1. DCs 
D1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_alpha", ident.2 = "DCs_untreated", sep = "_")
D1AlphaResponse$gene <- rownames(D1AlphaResponse)
D1AlphaResponse <- filter(D1AlphaResponse, p_val_adj < 0.05)
down <- filter(D1AlphaResponse, avg_log2FC < 0)
up <- filter(D1AlphaResponse, avg_log2FC >= 0)
saveRDS(D1AlphaResponse, "honours/results/FinalIndex/DEAnalysis/DendriticAlphaResponse.rds")

D1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_lambda", ident.2 = "DCs_untreated", sep = "_") 
D1LambdaResponse$gene <- rownames(D1LambdaResponse) 
D1LambdaResponse <- filter(D1LambdaResponse, p_val_adj < 0.05)
down <- filter(D1LambdaResponse, avg_log2FC < 0)
up <- filter(D1LambdaResponse, avg_log2FC >= 0)
saveRDS(D1LambdaResponse, "honours/results/FinalIndex/DEAnalysis/DendriticLambdaResponse.rds")

#2. mDCs
D2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "mDCs_alpha", ident.2 = "mDCs_untreated", sep = "_")
D2AlphaResponse$gene <- rownames(D2AlphaResponse)
D2AlphaResponse <- filter(D2AlphaResponse, p_val_adj < 0.05)
down <- filter(D2AlphaResponse, avg_log2FC < 0)
up <- filter(D2AlphaResponse, avg_log2FC >= 0)
saveRDS(D2AlphaResponse, "honours/results/FinalIndex/DEAnalysis/mdcsAlphaResponse.rds")

D2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "mDCs_lambda", ident.2 = "mDCs_untreated", sep = "_")
D2LambdaResponse$gene <- rownames(D2LambdaResponse) 
D2LambdaResponse <- filter(D2LambdaResponse, p_val_adj < 0.05)
down <- filter(D2LambdaResponse, avg_log2FC < 0)
up <- filter(D2LambdaResponse, avg_log2FC >= 0)
D2LambdaResponse <- data.frame(Gene = D2LambdaResponse$gene, Log2FoldChange = D2LambdaResponse$avg_log2FC)  
saveRDS(D2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/mdcsLambdaResponse.rds")

# 3. pDCs : only in alpha dataset, therefore no DEG
# D3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "pDCs_alpha", ident.2 = "pDCs_untreated", sep = "_")
# D3AlphaResponse$gene <- rownames(D3AlphaResponse)
# D3AlphaResponse <- filter(D3AlphaResponse, p_val_adj < 0.05)
# down <- filter(D3AlphaResponse, avg_log2FC < 0)
# up <- filter(D3AlphaResponse, avg_log2FC >= 0)
# D3AlphaResponse <- data.frame(Gene = D3AlphaResponse$gene, Log2FoldChange = D3AlphaResponse$avg_log2FC)
# saveRDS(D3AlphaResponse, "honours/results/FinalIndex/DEAnalysis/pdcsAlphaResponse.rds")
# 
# D3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "pDCs_lambda", ident.2 = "pDCs_untreated", sep = "_")
# D3LambdaResponse$gene <- rownames(D3LambdaResponse)
# D3LambdaResponse <- filter(D3LambdaResponse, p_val_adj < 0.05)
# down <- filter(D3LambdaResponse, avg_log2FC < 0)
# up <- filter(D3LambdaResponse, avg_log2FC >= 0)
# D3LambdaResponse <- data.frame(Gene = D3LambdaResponse$gene, Log2FoldChange = D3LambdaResponse$avg_log2FC)
# saveRDS(D3LambdaResponse, "honours/results/FinalIndex/DEAnalysis/pdcsLambdaResponse.csv")

# PLATELETS : 
PAlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_alpha", ident.2 = "platelets_untreated", sep = "_")
PAlphaResponse$gene <- rownames(PAlphaResponse)
PAlphaResponse <- filter(PAlphaResponse, p_val_adj < 0.05)
down <- filter(PAlphaResponse, avg_log2FC < 0)
up <- filter(PAlphaResponse, avg_log2FC >= 0)
saveRDS(PAlphaResponse, "honours/results/FinalIndex/DEAnalysis/plateletsAlphaResponse.rds")

PLambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_lambda", ident.2 = "platelets_untreated", sep = "_")
PLambdaResponse$gene <- rownames(PLambdaResponse) 
PLambdaResponse <- filter(PLambdaResponse, p_val_adj < 0.05)
down <- filter(PLambdaResponse, avg_log2FC < 0)
up <- filter(PLambdaResponse, avg_log2FC >= 0)
saveRDS(PLambdaResponse, "honours/results/FinalIndex/DEAnalysis/plateletsLambdaResponse.rds")

# T CELL TYPES : 
# 1. T cells
T1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_alpha", ident.2 = "T_untreated", sep = "_") 
T1AlphaResponse$gene <- rownames(T1AlphaResponse)
T1AlphaResponse <- filter(T1AlphaResponse, p_val_adj < 0.05)
down <- filter(T1AlphaResponse, avg_log2FC < 0)
up <- filter(T1AlphaResponse, avg_log2FC >= 0)
saveRDS(T1AlphaResponse, "honours/results/FinalIndex/DEAnalysis/TAlphaResponse.rds")

T1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_lambda", ident.2 = "T_untreated", sep = "_") 
T1LambdaResponse$gene <- rownames(T1LambdaResponse) 
T1LambdaResponse <- filter(T1LambdaResponse, p_val_adj < 0.05)
down <- filter(T1LambdaResponse, avg_log2FC < 0)
up <- filter(T1LambdaResponse, avg_log2FC >= 0)
saveRDS(T1LambdaResponse, "honours/results/FinalIndex/DEAnalysis/TLambdaResponse.rds")

# 2. CD4+ helper T 
T2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_helper_alpha", ident.2 = "CD4+_helper_untreated", sep = "_")  
T2AlphaResponse$gene <- rownames(T2AlphaResponse) 
T2AlphaResponse  <- filter(T2AlphaResponse , p_val_adj < 0.05)
down <- filter(T2AlphaResponse, avg_log2FC < 0)
up <- filter(T2AlphaResponse, avg_log2FC >= 0)
saveRDS(T2AlphaResponse, "honours/results/FinalIndex/DEAnalysis/CD4hAlphaResponse.rds")

T2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_helper_lambda", ident.2 = "CD4+_helper_untreated", sep = "_")  
T2LambdaResponse$gene <- rownames(T2LambdaResponse)
T2LambdaResponse <- filter(T2LambdaResponse, p_val_adj < 0.05)
down <- filter(T2LambdaResponse, avg_log2FC < 0)
up <- filter(T2LambdaResponse, avg_log2FC >= 0)
saveRDS(T2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/CD4hLambdaResponse.rds")

# 3. Naive CD8 T cells 
T3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "naive_CD8+ T_alpha", ident.2 = "naive_CD8+ T_untreated", sep = "_") 
T3AlphaResponse$gene <- rownames(T3AlphaResponse) 
T3AlphaResponse <- filter(T3AlphaResponse, p_val_adj < 0.05)
down <- filter(T3AlphaResponse, avg_log2FC < 0)
up <- filter(T3AlphaResponse, avg_log2FC >= 0)
saveRDS(T3AlphaResponse, "honours/results/FinalIndex/DEAnalysis/naivecd8AlphaResponse.rds")

T3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "naive_CD8+ T_lambda", ident.2 = "naive_CD8+ T_untreated", sep = "_")
T3LambdaResponse$gene <- rownames(T3LambdaResponse)
T3LambdaResponse <- filter(T3LambdaResponse, p_val_adj < 0.05)
down <- filter(T3LambdaResponse, avg_log2FC < 0)
up <- filter(T3LambdaResponse, avg_log2FC >= 0)
saveRDS(T3LambdaResponse, "honours/results/FinalIndex/DEAnalysis/naivecd8LambdaResponse.rds")
 
# 4. NKT 
T4AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NKT_alpha", ident.2 = "NKT_untreated", sep = "_")  
T4AlphaResponse$gene <- rownames(T4AlphaResponse) 
T4AlphaResponse <- filter(T4AlphaResponse, p_val_adj < 0.05)
down <- filter(T4AlphaResponse, avg_log2FC < 0)
up <- filter(T4AlphaResponse, avg_log2FC >= 0)
saveRDS(T4AlphaResponse, "honours/results/FinalIndex/DEAnalysis/NKTAlphaResponse.rds")

T4LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NKT_lambda", ident.2 = "NKT_untreated", sep = "_")
T4LambdaResponse$gene <- rownames(T4LambdaResponse)
T4LambdaResponse <- filter(T4LambdaResponse, p_val_adj < 0.05)
down <- filter(T4LambdaResponse, avg_log2FC < 0)
up <- filter(T4LambdaResponse, avg_log2FC >= 0)
saveRDS(T4LambdaResponse, "honours/results/FinalIndex/DEAnalysis/NKTLambdaResponse.rds")

# 5. cytotoxic cd8 t cells : 
T5AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cytotoxic_CD8+_T_alpha", ident.2 = "cytotoxic_CD8+_T_untreated", sep = "_")
T5AlphaResponse$gene <- rownames(T5AlphaResponse) 
T5AlphaResponse <- filter(T5AlphaResponse, p_val_adj < 0.05)
down <- filter(T5AlphaResponse, avg_log2FC < 0)
up <- filter(T5AlphaResponse, avg_log2FC >= 0)
saveRDS(T5AlphaResponse, "honours/results/FinalIndex/DEAnalysis/cytotoxiccd8AlphaResponse.rds")

T5LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cytotoxic_CD8+_T_lambda", ident.2 = "cytotoxic_CD8+_T_untreated", sep = "_") 
T5LambdaResponse$gene <- rownames(T5LambdaResponse)
T5LambdaResponse <- filter(T5LambdaResponse, p_val_adj < 0.05)
down <- filter(T5LambdaResponse, avg_log2FC < 0)
up <- filter(T5LambdaResponse, avg_log2FC >= 0)
saveRDS(T5LambdaResponse, "honours/results/FinalIndex/DEAnalysis/cytotoxiccd8LambdaResponse.rds")

# 6. Tregs : 
T6AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_alpha", ident.2 = "Tregs_untreated", sep = "_") 
T6AlphaResponse$gene <- rownames(T6AlphaResponse) 
T6AlphaResponse <- filter(T6AlphaResponse, p_val_adj < 0.05)
down <- filter(T6AlphaResponse, avg_log2FC < 0)
up <- filter(T6AlphaResponse, avg_log2FC >= 0)
saveRDS(T6AlphaResponse, "honours/results/FinalIndex/DEAnalysis/tregsAlphaResponse.rds")

T6LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_lambda", ident.2 = "Tregs_untreated", sep = "_") 
T6LambdaResponse$gene <- rownames(T6LambdaResponse)
T6LambdaResponse <- filter(T6LambdaResponse, p_val_adj < 0.05)
down <- filter(T6LambdaResponse, avg_log2FC < 0)
up <- filter(T6LambdaResponse, avg_log2FC >= 0)
saveRDS(T6LambdaResponse, "honours/results/FinalIndex/DEAnalysis/tregsLambdaResponse.rds")
 
# 7. CD4+ T cells : 
T7AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_T_alpha", ident.2 = "CD4+_T_untreated", sep = "_") 
T7AlphaResponse$gene <- rownames(T7AlphaResponse) 
T7AlphaResponse <- filter(T7AlphaResponse, p_val_adj < 0.05)
down <- filter(T7AlphaResponse, avg_log2FC < 0)
up <- filter(T7AlphaResponse, avg_log2FC >= 0)
saveRDS(T7AlphaResponse, "honours/results/FinalIndex/DEAnalysis/cd4AlphaResponse.rds")

T7LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_T_lambda", ident.2 = "CD4+_T_untreated", sep = "_")
T7LambdaResponse$gene <- rownames(T7LambdaResponse)
T7LambdaResponse <- filter(T7LambdaResponse, p_val_adj < 0.05)
down <- filter(T7LambdaResponse, avg_log2FC < 0)
up <- filter(T7LambdaResponse, avg_log2FC >= 0)
saveRDS(T7LambdaResponse, "honours/results/FinalIndex/DEAnalysis/cd4LambdaResponse.rds")

# 8. NK cells 
T8AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_alpha", ident.2 = "NK_untreated", sep = "_")
T8AlphaResponse$gene <- rownames(T8AlphaResponse) 
T8AlphaResponse <- filter(T8AlphaResponse, p_val_adj < 0.05)
down <- filter(T8AlphaResponse, avg_log2FC < 0)
up <- filter(T8AlphaResponse, avg_log2FC >= 0)
saveRDS(T8AlphaResponse, "honours/results/FinalIndex/DEAnalysis/NKAlphaResponse.rds")

T8LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_lambda", ident.2 = "NK_untreated", sep = "_") 
T8LambdaResponse$gene <- rownames(T8LambdaResponse)
T8LambdaResponse <- filter(T8LambdaResponse, p_val_adj < 0.05)
down <- filter(T8LambdaResponse, avg_log2FC < 0)
up <- filter(T8LambdaResponse, avg_log2FC >= 0)
saveRDS(T8LambdaResponse, "honours/results/FinalIndex/DEAnalysis/NKLambdaResponse.rds")

# B CELLS 
BAlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_alpha", ident.2 = "B_untreated", sep = "_") 
BAlphaResponse$gene <- rownames(BAlphaResponse) 
BAlphaResponse <- filter(BAlphaResponse, p_val_adj < 0.05)
down <- filter(BAlphaResponse, avg_log2FC < 0)
up <- filter(BAlphaResponse, avg_log2FC >= 0)
saveRDS(BAlphaResponse, "honours/results/FinalIndex/DEAnalysis/BAlphaResponse.rds")

BLambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_lambda", ident.2 = "B_untreated", sep = "_")
BLambdaResponse$gene <- rownames(BLambdaResponse)
BLambdaResponse <- filter(BLambdaResponse, p_val_adj < 0.05)
down <- filter(BLambdaResponse, avg_log2FC < 0)
up <- filter(BLambdaResponse, avg_log2FC >= 0)
saveRDS(BLambdaResponse, "honours/results/FinalIndex/DEAnalysis/BLambdaResponse.rds")



##### [2.2] DE with FindMarkers() on broad cell types () #####
treatment <- readRDS("honours/results/FinalIndex/adjtreatment.rds")
DefaultAssay(treatment) <- "RNA"

# change annotations to create larger groups to test DE on : 
TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'myeloid',
                                   '1' = 'T',
                                   '2' = 'myeloid',
                                   '3' = 'T',
                                   '4' = 'myeloid',
                                   '5' = 'T',
                                   '6' = 'B',
                                   '7' = 'T',
                                   '8' = 'DCs',
                                   '9' = 'T',
                                   '10' = 'T',
                                   '11' = 'T',
                                   '12' = 'platelets',
                                   '13' = 'unknown',
                                   '14' = 'DCs', 
                                   '15' = 'T',
                                   '16' = 'T',
                                   '17' = 'DCs')


TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$treatment, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the idents to that column 
# Now we can easily refer to a cell & treatment type, to see the options : 
levels(TreatmentAnnotated) # list the options to perform DE on :

# MYELOID CELL TYPES : 
system.time(MAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "myeloid_alpha", ident.2 = "myeloid_untreated", sep = "_"))
MAResponse$gene <- rownames(MAResponse) # puts gene names into a column
MAResponse <- filter(MAResponse, p_val_adj < 0.05)
down <- filter(MAResponse, avg_log2FC < 0)
up <- filter(MAResponse, avg_log2FC >= 0)
saveRDS(MAResponse, "honours/results/FinalIndex/DEAnalysis/MyeloidAlphaResponse.rds")

system.time(MLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "myeloid_lambda", ident.2 = "myeloid_untreated", sep = "_"))
MLResponse$gene <- rownames(MLResponse) # puts gene names into a column
MLResponse <- filter(MLResponse, p_val_adj < 0.05)
down <- filter(MLResponse, avg_log2FC < 0)
up <- filter(MLResponse, avg_log2FC >= 0)
saveRDS(MLResponse, "honours/results/FinalIndex/DEAnalysis/MyeloidLambdaResponse.rds")

# Dendritic cell types : 
system.time(DAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_alpha", ident.2 = "DCs_untreated", sep = "_"))
DAResponse$gene <- rownames(DAResponse) # puts gene names into a column
DAResponse <- filter(DAResponse, p_val_adj < 0.05)
down <- filter(DAResponse, avg_log2FC < 0)
up <- filter(DAResponse, avg_log2FC >= 0)
saveRDS(DAResponse, "honours/results/FinalIndex/DEAnalysis/DendriticAlphaResponse.rds")

system.time(DLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_lambda", ident.2 = "DCs_untreated", sep = "_"))
DLResponse$gene <- rownames(DLResponse) # puts gene names into a column
DLResponse <- filter(DLResponse, p_val_adj < 0.05)
down <- filter(DLResponse, avg_log2FC < 0)
up <- filter(DLResponse, avg_log2FC >= 0)
saveRDS(DLResponse, "honours/results/FinalIndex/DEAnalysis/DendriticLambdaResponse.rds")


# T cells : 
system.time(TAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_alpha", ident.2 = "T_untreated", sep = "_"))
TAResponse$gene <- rownames(TAResponse) # puts gene names into a column
TAResponse <- filter(TAResponse, p_val_adj < 0.05)
down <- filter(TAResponse, avg_log2FC < 0)
up <- filter(TAResponse, avg_log2FC >= 0)
saveRDS(TAResponse, "honours/results/FinalIndex/DEAnalysis/allTAlphaResponse.rds")

system.time(TLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_lambda", ident.2 = "T_untreated", sep = "_"))
TLResponse$gene <- rownames(TLResponse) # puts gene names into a column
TLResponse <- filter(TLResponse, p_val_adj < 0.05)
down <- filter(TLResponse, avg_log2FC < 0)
up <- filter(TLResponse, avg_log2FC >= 0)
saveRDS(TLResponse, "honours/results/FinalIndex/DEAnalysis/allTLambdaResponse.rds")


##### [2.3] Looking for common & unique DEGs #####
lambda_datasets <- list(M1LambdaResponse, M2LambdaResponse, T1LambdaResponse, T2LambdaResponse, 
  T3LambdaResponse, T4LambdaResponse, T5LambdaResponse, T6LambdaResponse, 
  T8LambdaResponse, BLambdaResponse)

alpha_datasets <- list(M1AlphaResponse, M2AlphaResponse, T1AlphaResponse, T2AlphaResponse, 
  T3AlphaResponse, T4AlphaResponse, T5AlphaResponse, T6AlphaResponse, 
  T8AlphaResponse, BAlphaResponse)

# Initialize a list to store the number of unique DEGs for each comparison
unique_degs_counts <- list()

# Loop through each lambda dataset
for (i in 1:length(lambda_datasets)) {            # Extract gene names from the current lambda dataset
  lambda_genes <- lambda_datasets[[i]]$gene       # Extract gene names from the corresponding alpha dataset
  alpha_genes <- alpha_datasets[[i]]$gene         # Find unique DEGs in the lambda dataset but not in the alpha dataset
  unique_lambda_genes <- setdiff(lambda_genes, alpha_genes) # Get the count of unique DEGs for this comparison
  num_unique_lambda_genes <- length(unique_lambda_genes)    # Store the count in the list
  unique_degs_counts[[i]] <- num_unique_lambda_genes        # Print the number of unique DEGs for this comparison
  cat("Number of unique DEGs in lambda dataset", i, "over alpha dataset", i, ":", num_unique_lambda_genes, "\n")
}



# [c] Store information in excel file 
B_cells <- list('AlphaResponse' = B_cells_alpha, 'LambdaResponse' = B_cells_lambda,'common' = B_cells_common_dataframe)
openxlsx::write.xlsx(B_cells, file = "honours/results/FinalIndex/DEAnalysis/B_cells.xlsx")

##### [3] ClusterProfiler GO analysis #####

# 1 : install and load annotation for desired organism (human) : org.Hs.eg.db
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# 2 : have to convert DE gene list to contain entrezgene_ID BIOMART
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Create a Mart object for the human genome
# listDatasets(mart) # List available datasets for the Mart (214!)

# 3:  Retrieve gene information for human genes :
gene_ID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
                 mart = mart) 

# # 4: create a background universe from the 2000 most variable genes 
# # Assuming 'mySeurat' is your integrated Seurat object
# # Access the integrated assay slot
# integrated_assay <- TreatmentAnnotated@assays[["integrated"]]
# rna_assay <- rownames(TreatmentAnnotated@assays[['RNA']])
# # Get the 2000 most variable genes
# top_variable_genes <- rownames(integrated_assay)[order(integrated_assay, decreasing = TRUE)[1:2000]]
# universe <- select(org.Hs.eg.db, keys = rna_assay, keytype = "SYMBOL", columns = "ENTREZID")
# universe <- universe$ENTREZID


#### REPITITION STARTS HERE : 
# 4 : merge this with DEGs datasets : 
# myeloid 
M1AlphaResponse <- unique(merge(M1AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name"))
M1LambdaResponse <- unique(merge(M1LambdaResponse, gene_ID[,c(2,3)],by.x = "gene", by.y = "external_gene_name"))
M2AlphaResponse <- merge(M2AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
M2LambdaResponse <- merge(M2LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
MAResponse <- merge(MAResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
MLResponse <- merge(MLResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

# dendritic 
D1AlphaResponse <- merge(D1AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
D1LambdaResponse <- merge(D1LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
D2AlphaResponse <- merge(D2AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
D2LambdaResponse <- merge(D2LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
# D3AlphaResponse <- merge(D3AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
# D3LambdaResponse <- merge(D3LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
DAResponse <- merge(DAResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
DLResponse <- merge(DLResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

# platelets : 
PAlphaResponse <- merge(PAlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
PLambdaResponse <- merge(PLambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

# T cells : 
T1AlphaResponse <- merge(T1AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T1LambdaResponse <- merge(T1LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T2AlphaResponse <- merge(T2AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T2LambdaResponse <- merge(T2LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T3AlphaResponse <- merge(T3AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T3LambdaResponse <- merge(T3LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T4AlphaResponse <- merge(T4AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T4LambdaResponse <- merge(T4LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T5AlphaResponse <- merge(T5AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T5LambdaResponse <- merge(T5LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T6AlphaResponse <- merge(T6AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T6LambdaResponse <- merge(T6LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T7AlphaResponse <- merge(T7AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T7LambdaResponse <- merge(T7LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T8AlphaResponse <- merge(T8AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
T8LambdaResponse <- merge(T8LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
TAResponse <- merge(TAResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
TLResponse <- merge(TLResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")


# B cells 
BAlphaResponse <- merge(BAlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
BLambdaResponse <- unique(merge(BLambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name"))

# 5 : create object to perform GO analysis on, just the entrezgene  id column : 
# myeloid 
M1AlphaDE <- na.omit(M1AlphaResponse$entrezgene_id)
M1LambdaDE <- na.omit(M1LambdaResponse$entrezgene_id)
M2AlphaDE <- na.omit(M2AlphaResponse$entrezgene_id)
M2LambdaDE <- na.omit(M2LambdaResponse$entrezgene_id)
MAlphaDE <- na.omit(MAResponse$entrezgene_id)
MLambdaDE <- na.omit(MLResponse$entrezgene_id)

# dendritic 
D1AlphaDE <- na.omit(D1AlphaResponse$entrezgene_id)
D1LambdaDE <- na.omit(D1LambdaResponse$entrezgene_id)
D2AlphaDE <- na.omit(D2AlphaResponse$entrezgene_id)
D2LambdaDE <- na.omit(D2LambdaResponse$entrezgene_id)
# D3AlphaDE <- D3AlphaResponse$entrezgene_id
# D3LambdaDE <- D3LambdaResponse$entrezgene_id
DAlphaDE <- na.omit(DAResponse$entrezgene_id)
DLambdaDE <- na.omit(DLResponse$entrezgene_id)

# platelets 
PAlphaDE <- na.omit(PAlphaResponse$entrezgene_id)
PLambdaDE <- na.omit(PLambdaResponse$entrezgene_id)

# T cells 
T1AlphaDE <- na.omit(T1AlphaResponse$entrezgene_id)
T1LambdaDE <- na.omit(T1LambdaResponse$entrezgene_id)
T2AlphaDE <- na.omit(T2AlphaResponse$entrezgene_id)
T2LambdaDE <- na.omit(T2LambdaResponse$entrezgene_id)
T3AlphaDE <- na.omit(T3AlphaResponse$entrezgene_id)
T3LambdaDE <- na.omit(T3LambdaResponse$entrezgene_id)
T4AlphaDE <- na.omit(T4AlphaResponse$entrezgene_id)
T4LambdaDE <- na.omit(T4LambdaResponse$entrezgene_id)
T5AlphaDE <- na.omit(T5AlphaResponse$entrezgene_id)
T5LambdaDE <- na.omit(T5LambdaResponse$entrezgene_id)
T6AlphaDE <- na.omit(T6AlphaResponse$entrezgene_id)
T6LambdaDE <- na.omit(T6LambdaResponse$entrezgene_id)
T7AlphaDE <- na.omit(T7AlphaResponse$entrezgene_id)
T7LambdaDE <- na.omit(T7LambdaResponse$entrezgene_id)
T8AlphaDE <- na.omit(T8AlphaResponse$entrezgene_id)
T8LambdaDE <- na.omit(T8LambdaResponse$entrezgene_id)
TAlphaDE <- na.omit(TAResponse$entrezgene_id)
TLambdaDE <- na.omit(TLResponse$entrezgene_id)

# B cells 
BAlphaDE <- na.omit(BAlphaResponse$entrezgene_id)
BLambdaDE <- na.omit(BLambdaResponse$entrezgene_id)

# 6 : perform GO analysis 
?enrichGO



# Myeloid  
M1AlphaEGO <- enrichGO(gene = M1AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.05)
M1AlphaEGO <- filter(M1AlphaEGO, M1AlphaEGO@result$p.adjust < 0.05)
M1LambdaEGO <- enrichGO(gene = M1LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
M1LambdaEGO <- filter(M1LambdaEGO, M1LambdaEGO@result$p.adjust < 0.05)
M2AlphaEGO <- enrichGO(gene = M2AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
M2AlphaEGO <- filter(M2AlphaEGO, M2AlphaEGO@result$p.adjust < 0.05)
M2LambdaEGO <- enrichGO(gene = M2LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
M2LambdaEGO <- filter(M2LambdaEGO, M2LambdaEGO@result$p.adjust < 0.05)
MAlphaEGO <- enrichGO(gene = MAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
MAlphaEGO <- filter(MAlphaEGO, MAlphaEGO@result$p.adjust < 0.05)
MLambdaEGO <- enrichGO(gene = MLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
MLambdaEGO <- filter(MLambdaEGO, MLambdaEGO@result$p.adjust < 0.05)

# Dendritic cells ; 
D1AlphaEGO <- enrichGO(gene = D1AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
D1AlphaEGO <- filter(D1AlphaEGO, D1AlphaEGO@result$p.adjust < 0.05)
# D1LambdaEGO <- enrichGO(gene = D1LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D1LambdaEGO <- filter(D1LambdaEGO, D1LambdaEGO@result$p.adjust < 0.05)
D2AlphaEGO <- enrichGO(gene = D2AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
D2AlphaEGO <- filter(D2AlphaEGO, D2AlphaEGO@result$p.adjust < 0.05)
# D2LambdaEGO <- enrichGO(gene = D2LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D2LambdaEGO <- filter(D2LambdaEGO, D2LambdaEGO@result$p.adjust < 0.05)
# D3AlphaEGO <- enrichGO(gene = D3AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D3AlphaEGO <- filter(D3AlphaEGO, D3AlphaEGO@result$p.adjust < 0.05)
# D3LambdaEGO <- enrichGO(gene = D3LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D3LambdaEGO <- filter(D3LambdaEGO, D3LambdaEGO@result$p.adjust < 0.05)
DAlphaEGO <- enrichGO(gene = DAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
DAlphaEGO <- filter(DAlphaEGO, DAlphaEGO@result$p.adjust < 0.05)


# Platelets 
PAlphaEGO <- enrichGO(gene = PAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
PAlphaEGO <- filter(PAlphaEGO, PAlphaEGO@result$p.adjust < 0.05)
# PLambdaEGO <- enrichGO(gene = PLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# PLambdaEGO <- filter(PLambdaEGO, PLambdaEGO@result$p.adjust < 0.05)

# T cells  
T1AlphaEGO <- enrichGO(gene = T1AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T1AlphaEGO <- filter(T1AlphaEGO, T1AlphaEGO@result$p.adjust < 0.05)
T1LambdaEGO <- enrichGO(gene = T1LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T1LambdaEGO <- filter(T1LambdaEGO, T1LambdaEGO@result$p.adjust < 0.05)
T2AlphaEGO <- enrichGO(gene = T2AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T2AlphaEGO <- filter(T2AlphaEGO, T2AlphaEGO@result$p.adjust < 0.05)
T2LambdaEGO <- enrichGO(gene = T2LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T2LambdaEGO <- filter(T2LambdaEGO, T2LambdaEGO@result$p.adjust < 0.05)
T3AlphaEGO <- enrichGO(gene = T3AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T3AlphaEGO <- filter(T3AlphaEGO, T3AlphaEGO@result$p.adjust < 0.05)
T3LambdaEGO <- enrichGO(gene = T3LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T3LambdaEGO <- filter(T3LambdaEGO, T3LambdaEGO@result$p.adjust < 0.05)
T4AlphaEGO <- enrichGO(gene = T4AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T4AlphaEGO <- filter(T4AlphaEGO, T4AlphaEGO@result$p.adjust < 0.05)
T4LambdaEGO <- enrichGO(gene = T4LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T4LambdaEGO <- filter(T4LambdaEGO, T4LambdaEGO@result$p.adjust < 0.05)
T5AlphaEGO <- enrichGO(gene = T5AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T5AlphaEGO <- filter(T5AlphaEGO, T5AlphaEGO@result$p.adjust < 0.05)
T5LambdaEGO <- enrichGO(gene = T5LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T5LambdaEGO <- filter(T5LambdaEGO, T5LambdaEGO@result$p.adjust < 0.05)
T6AlphaEGO <- enrichGO(gene = T6AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T6AlphaEGO <- filter(T6AlphaEGO, T6AlphaEGO@result$p.adjust < 0.05)
T6LambdaEGO <- enrichGO(gene = T6LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T6LambdaEGO <- filter(T6LambdaEGO, T6LambdaEGO@result$p.adjust < 0.05)
T7AlphaEGO <- enrichGO(gene = T7AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T7AlphaEGO <- filter(T7AlphaEGO, T7AlphaEGO@result$p.adjust < 0.05)
# T7LambdaEGO <- enrichGO(gene = T7LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# T7LambdaEGO <- filter(T7LambdaEGO, T7LambdaEGO@result$p.adjust < 0.05)
T8AlphaEGO <- enrichGO(gene = T8AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T8AlphaEGO <- filter(T8AlphaEGO, T8AlphaEGO@result$p.adjust < 0.05)
T8LambdaEGO <- enrichGO(gene = T8LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T8LambdaEGO <- filter(T8LambdaEGO, T8LambdaEGO@result$p.adjust < 0.05)
TAlphaEGO <- enrichGO(gene = TAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
TAlphaEGO <- filter(TAlphaEGO, TAlphaEGO@result$p.adjust < 0.05)
TLambdaEGO <- enrichGO(gene = TLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
TLambdaEGO <- filter(TLambdaEGO, TLambdaEGO@result$p.adjust < 0.05)

# B CELLS 
BAlphaEGO <- enrichGO(gene = BAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.05)
BAlphaEGO <- filter(BAlphaEGO, BAlphaEGO@result$p.adjust < 0.05)
BLambdaEGO <- enrichGO(gene = BLambdaDE, OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",  pvalueCutoff = 0.05)
BLambdaEGO <- filter(BLambdaEGO, BLambdaEGO@result$qvalue < 0.10)

# save enrich GO output : 
saveRDS(M1AlphaEGO, "honours/results/FinalIndex/DEAnalysis/MonocytesAlphaEGO.rds")
saveRDS(M1LambdaEGO, "honours/results/FinalIndex/DEAnalysis/MonocytesLambdaEGO.rds")
saveRDS(M2AlphaEGO, "honours/results/FinalIndex/DEAnalysis/NeutrophilsAlphaEGO.rds")
saveRDS(M2LambdaEGO, "honours/results/FinalIndex/DEAnalysis/NeutrophilsLambdaEGO.rds")
saveRDS(MAlphaEGO, "honours/results/FinalIndex/DEAnalysis/MyeloidAlphaEGO.rds")
saveRDS(MLambdaEGO, "honours/results/FinalIndex/DEAnalysis/MyeloidLambdaEGO.rds")

saveRDS(D1AlphaEGO, "honours/results/FinalIndex/DEAnalysis/DendriticAlphaEGO.rds")
saveRDS(D2AlphaEGO, "honours/results/FinalIndex/DEAnalysis/mDendriticAlphaEGO.rds")
saveRDS(DAlphaEGO, "honours/results/FinalIndex/DEAnalysis/AllDendriticAlphaEGO.rds")

saveRDS(PAlphaEGO, "honours/results/FinalIndex/DEAnalysis/PlateletsAlphaEGO.rds")

saveRDS(T1AlphaEGO, "honours/results/FinalIndex/DEAnalysis/TAlphaEGO.rds")
saveRDS(T1LambdaEGO, "honours/results/FinalIndex/DEAnalysis/TLambdaEGO.rds")
saveRDS(T2AlphaEGO, "honours/results/FinalIndex/DEAnalysis/CD4hAlphaEGO.rds")
saveRDS(T2LambdaEGO, "honours/results/FinalIndex/DEAnalysis/CD4hLambdaEGO.rds")
saveRDS(T3AlphaEGO, "honours/results/FinalIndex/DEAnalysis/naiveCD8AlphaEGO.rds")
saveRDS(T3LambdaEGO, "honours/results/FinalIndex/DEAnalysis/naiveCD8LambdaEGO.rds")
saveRDS(T4AlphaEGO, "honours/results/FinalIndex/DEAnalysis/NKTAlphaEGO.rds")
saveRDS(T4LambdaEGO, "honours/results/FinalIndex/DEAnalysis/NKTLambdaEGO.rds")
saveRDS(T5AlphaEGO, "honours/results/FinalIndex/DEAnalysis/cytoCD8AlphaEGO.rds")
saveRDS(T5LambdaEGO, "honours/results/FinalIndex/DEAnalysis/cytoCD8LambdaEGO.rds")
saveRDS(T6AlphaEGO, "honours/results/FinalIndex/DEAnalysis/TregslphaEGO.rds")
saveRDS(T6LambdaEGO, "honours/results/FinalIndex/DEAnalysis/TregsLambdaEGO.rds")
saveRDS(T7AlphaEGO, "honours/results/FinalIndex/DEAnalysis/CD4AlphaEGO.rds")
#saveRDS(T7LambdaEGO, "honours/results/FinalIndex/DEAnalysis/CD4LambdaEGO.rds")
saveRDS(T8AlphaEGO, "honours/results/FinalIndex/DEAnalysis/NKAlphaEGO.rds")
saveRDS(T8LambdaEGO, "honours/results/FinalIndex/DEAnalysis/NKLambdaEGO.rds")
saveRDS(TAlphaEGO, "honours/results/FinalIndex/DEAnalysis/AllTAlphaEGO.rds")
saveRDS(TLambdaEGO, "honours/results/FinalIndex/DEAnalysis/AllTLambdaEGO.rds")

saveRDS(BAlphaEGO, "honours/results/FinalIndex/DEAnalysis/BAlphaEGO.rds")
saveRDS(BLambdaEGO, "honours/results/FinalIndex/DEAnalysis/BLambdaEGO.rds")

# Then, proceed with enrichGO
# 7 : of the descriptions found in alpha & lambda set, subset according to unique to alpha and lambda and the intersection 
# Myeloid : 
M1_common <- intersect(M1AlphaEGO@result$Description, M1LambdaEGO@result$Description )                # identify common genes
M1_alpha <- M1AlphaEGO@result$Description[!(M1AlphaEGO@result$Description %in% M1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
M1_lambda <- M1LambdaEGO@result$Description[!(M1LambdaEGO@result$Description %in% M1AlphaEGO@result$Description)]
M1_lambda <- M1LambdaEGO[M1LambdaEGO@result$Description %in% M1_lambda, ] 
M1_lambda_list <- M1_lambda[, c("ID", "p.adjust")] 
write_tsv(M1_lambda_list, "honours/results/FinalIndex/GOAnalysis/M1_lambda_list.tsv")

M2_common <- intersect(M2AlphaEGO@result$Description, M2LambdaEGO@result$Description )                # identify common genes
M2_alpha <- M2AlphaEGO@result$Description[!(M2AlphaEGO@result$Description %in% M2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
M2_lambda <- M2LambdaEGO@result$Description[!(M2LambdaEGO@result$Description %in% M2AlphaEGO@result$Description)]
M2_lambda <- M2LambdaEGO[M2LambdaEGO@result$Description %in% M2_lambda, ] 
M2_lambda_list <- M2_lambda[, c("ID", "p.adjust")] 
write_tsv(M2_lambda_list, "honours/results/FinalIndex/GOAnalysis/M2_lambda_list.tsv")

M_common <- intersect(MAlphaEGO@result$Description, MLambdaEGO@result$Description )                # identify common genes
M_alpha <- MAlphaEGO@result$Description[!(MAlphaEGO@result$Description %in% MLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
M_lambda <- MLambdaEGO@result$Description[!(MLambdaEGO@result$Description %in% MAlphaEGO@result$Description)]
M_lambda <- MLambdaEGO[MLambdaEGO@result$Description %in% M_lambda, ] 
M_lambda_list <- M_lambda[, c("ID", "p.adjust")] 
write_tsv(M_lambda_list, "honours/results/FinalIndex/GOAnalysis/M_lambda_list.tsv")

# Dendritic cells : 
# D1_common <- intersect(D1AlphaEGO@result$Description, D1LambdaEGO@result$Description )                # identify common genes
# D1_alpha <- D1AlphaEGO@result$Description[!(D1AlphaEGO@result$Description %in% D1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# D1_lambda <- D1LambdaEGO@result$Description[!(D1LambdaEGO@result$Description %in% D1AlphaEGO@result$Description)]
# D1_lambda <- D1LambdaEGO[D1LambdaEGO@result$Description %in% D1_lambda, ] 
# D1_lambda_list <- D1_lambda[, c("ID", "p.adjust")] 
# write_tsv(D1_lambda_list, "honours/results/FinalIndex/GOAnalysis/D1_lambda_list.tsv")

# D2_common <- intersect(D2AlphaEGO@result$Description, D2LambdaEGO@result$Description )                # identify common genes
# D2_alpha <- D2AlphaEGO@result$Description[!(D2AlphaEGO@result$Description %in% D2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# D2_lambda <- D2LambdaEGO@result$Description[!(D2LambdaEGO@result$Description %in% D2AlphaEGO@result$Description)]
# D2_lambda <- D2LambdaEGO[D2LambdaEGO@result$Description %in% D2_lambda, ] 
# D2_lambda_list <- D2_lambda[, c("ID", "p.adjust")] 
# write_tsv(D2_lambda_list, "honours/results/FinalIndex/GOAnalysis/D2_lambda_list.tsv")

# D3_common <- intersect(D3AlphaEGO@result$Description, D3LambdaEGO@result$Description )                # identify common genes
# D3_alpha <- D3AlphaEGO@result$Description[!(D3AlphaEGO@result$Description %in% D3LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# D3_lambda <- D3LambdaEGO@result$Description[!(D3LambdaEGO@result$Description %in% D3AlphaEGO@result$Description)]
# D3_lambda <- D3LambdaEGO[D3LambdaEGO@result$Description %in% D3_lambda, ] 
# D3_lambda_list <- D3_lambda[, c("ID", "p.adjust")] 
# write_tsv(D3_lambda_list, "honours/results/FinalIndex/GOAnalysis/D3_lambda_list.tsv")

# D_common <- intersect(DAlphaEGO@result$Description, DLambdaEGO@result$Description )                # identify common genes
# D_alpha <- DAlphaEGO@result$Description[!(DAlphaEGO@result$Description %in% DLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# D_lambda <- DLambdaEGO@result$Description[!(DLambdaEGO@result$Description %in% DAlphaEGO@result$Description)]
# D_lambda <- DLambdaEGO[DLambdaEGO@result$Description %in% D_lambda, ] 
# D_lambda_list <- D_lambda[, c("ID", "p.adjust")] 
# write_tsv(D_lambda_list, "honours/results/FinalIndex/GOAnalysis/D_lambda_list.tsv")

# Platelets : 
# P_common <- intersect(PAlphaEGO@result$Description, PLambdaEGO@result$Description )                # identify common genes
# P_alpha <- PAlphaEGO@result$Description[!(PAlphaEGO@result$Description %in% PLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# P_lambda <- PLambdaEGO@result$Description[!(PLambdaEGO@result$Description %in% PAlphaEGO@result$Description)]
# P_lambda <- PLambdaEGO[PLambdaEGO@result$Description %in% P_lambda, ] 
# P_lambda_list <- P_lambda[, c("ID", "p.adjust")] 
# write_tsv(P_lambda_list, "honours/results/FinalIndex/GOAnalysis/P_lambda_list.tsv")

# T cells : 
T1_common <- intersect(T1AlphaEGO@result$Description, T1LambdaEGO@result$Description )                # identify common genes
T1_alpha <- T1AlphaEGO@result$Description[!(T1AlphaEGO@result$Description %in% T1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T1_lambda <- T1LambdaEGO@result$Description[!(T1LambdaEGO@result$Description %in% T1AlphaEGO@result$Description)]
T1_lambda <- T1LambdaEGO[T1LambdaEGO@result$Description %in% T1_lambda, ]
T1_lambda_list <- T1_lambda[, c("ID", "p.adjust")]
write_tsv(T1_lambda_list, "honours/results/FinalIndex/GOAnalysis/T1_lambda_list.tsv")

T2_common <- intersect(T2AlphaEGO@result$Description, T2LambdaEGO@result$Description )                # identify common genes
T2_alpha <- T2AlphaEGO@result$Description[!(T2AlphaEGO@result$Description %in% T2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T2_lambda <- T2LambdaEGO@result$Description[!(T2LambdaEGO@result$Description %in% T2AlphaEGO@result$Description)]
T2_lambda <- T2LambdaEGO[T2LambdaEGO@result$Description %in% T2_lambda, ]
T2_lambda_list <- T2_lambda[, c("ID", "p.adjust")]
write_tsv(T2_lambda_list, "honours/results/FinalIndex/GOAnalysis/T2_lambda_list.tsv")

T3_common <- intersect(T3AlphaEGO@result$Description, T3LambdaEGO@result$Description )                # identify common genes
T3_alpha <- T3AlphaEGO@result$Description[!(T3AlphaEGO@result$Description %in% T3LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T3_lambda <- T3LambdaEGO@result$Description[!(T3LambdaEGO@result$Description %in% T3AlphaEGO@result$Description)]
T3_lambda <- T3LambdaEGO[T3LambdaEGO@result$Description %in% T3_lambda, ]
T3_lambda_list <- T3_lambda[, c("ID", "p.adjust")]
write_tsv(T3_lambda_list, "honours/results/FinalIndex/GOAnalysis/T3_lambda_list.tsv")

T4_common <- intersect(T4AlphaEGO@result$Description, T4LambdaEGO@result$Description )                # identify common genes
T4_alpha <- T4AlphaEGO@result$Description[!(T4AlphaEGO@result$Description %in% T4LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T4_lambda <- T4LambdaEGO@result$Description[!(T4LambdaEGO@result$Description %in% T4AlphaEGO@result$Description)]
T4_lambda <- T4LambdaEGO[T4LambdaEGO@result$Description %in% T4_lambda, ] 
T4_lambda_list <- T4_lambda[, c("ID", "p.adjust")] 
write_tsv(T4_lambda_list, "honours/results/FinalIndex/GOAnalysis/T4_lambda_list.tsv")

T5_common <- intersect(T5AlphaEGO@result$Description, T5LambdaEGO@result$Description )                # identify common genes
T5_alpha <- T5AlphaEGO@result$Description[!(T5AlphaEGO@result$Description %in% T5LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T5_lambda <- T5LambdaEGO@result$Description[!(T5LambdaEGO@result$Description %in% T5AlphaEGO@result$Description)]
T5_lambda <- T5LambdaEGO[T5LambdaEGO@result$Description %in% T5_lambda, ] 
T5_lambda_list <- T5_lambda[, c("ID", "p.adjust")] 
write_tsv(T5_lambda_list, "honours/results/FinalIndex/GOAnalysis/T5_lambda_list.tsv")

T6_common <- intersect(T6AlphaEGO@result$Description, T6LambdaEGO@result$Description )                # identify common genes
T6_alpha <- T6AlphaEGO@result$Description[!(T6AlphaEGO@result$Description %in% T6LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T6_lambda <- T6LambdaEGO@result$Description[!(T6LambdaEGO@result$Description %in% T6AlphaEGO@result$Description)]
T6_lambda <- T6LambdaEGO[T6LambdaEGO@result$Description %in% T6_lambda, ] 
T6_lambda_list <- T6_lambda[, c("ID", "p.adjust")] 
write_tsv(T6_lambda_list, "honours/results/FinalIndex/GOAnalysis/T6_lambda_list.tsv")

# T7_common <- intersect(T7AlphaEGO@result$Description, T7LambdaEGO@result$Description )                # identify common genes
# T7_alpha <- T7AlphaEGO@result$Description[!(T7AlphaEGO@result$Description %in% T7LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# T7_lambda <- T7LambdaEGO@result$Description[!(T7LambdaEGO@result$Description %in% T7AlphaEGO@result$Description)]
# T7_lambda <- T7LambdaEGO[T7LambdaEGO@result$Description %in% T7_lambda, ] 
# T7_lambda_list <- T7_lambda[, c("ID", "p.adjust")] 
# write_tsv(T7_lambda_list, "honours/results/FinalIndex/GOAnalysis/T7_lambda_list.tsv")

T8_common <- intersect(T8AlphaEGO@result$Description, T8LambdaEGO@result$Description )                # identify common genes
T8_alpha <- T8AlphaEGO@result$Description[!(T8AlphaEGO@result$Description %in% T8LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T8_lambda <- T8LambdaEGO@result$Description[!(T8LambdaEGO@result$Description %in% T8AlphaEGO@result$Description)]
T8_lambda <- T8LambdaEGO[T8LambdaEGO@result$Description %in% T8_lambda, ] 
T8_lambda_list <- T8_lambda[, c("ID", "p.adjust")] 
write_tsv(T8_lambda_list, "honours/results/FinalIndex/GOAnalysis/T8_lambda_list.tsv")

T_common <- intersect(TAlphaEGO@result$Description, TLambdaEGO@result$Description )                # identify common genes
T_alpha <- TAlphaEGO@result$Description[!(TAlphaEGO@result$Description %in% TLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T_lambda <- TLambdaEGO@result$Description[!(TLambdaEGO@result$Description %in% TAlphaEGO@result$Description)]
T_lambda <- TLambdaEGO[TLambdaEGO@result$Description %in% T_lambda, ] 
T_lambda_list <- T_lambda[, c("ID", "p.adjust")] 
write_tsv(T_lambda_list, "honours/results/FinalIndex/GOAnalysis/T_lambda_list.tsv")


# BCELLS : 
B_common <- intersect(BAlphaEGO@result$Description, BLambdaEGO@result$Description )                # identify common genes
B_alpha <- BAlphaEGO@result$Description[!(BAlphaEGO@result$Description %in% BLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
B_lambda <- BLambdaEGO@result$Description[!(BLambdaEGO@result$Description %in% BAlphaEGO@result$Description)]
B_lambda <- BLambdaEGO[BLambdaEGO@result$Description %in% B_lambda, ] 
B_lambda_list <- B_lambda[, c("ID", "p.adjust")] 
write_tsv(B_lambda_list, "honours/results/FinalIndex/GOAnalysis/B_lambda_list.tsv")

# Saving the data to an excel sheet to annotate genes : 
GOmyeloid <- list(# Myeloid cells : 
                'MonocytesCommon' = M1_common, 'MonocytesAlpha' = M1_alpha,'MonocytesLambda' = M1_lambda,
                'NeutrophilsCommon' = M2_common, 'NeutrophilsAlpha' = M2_alpha,'NeutrophilsLambda' = M2_lambda, 
                'MyeloidCommon' = M_common, 'MyeloidAlpha' = M_alpha,'MyeloidLambda' = M_lambda)
# GOdendritic <- list(# Dendritic cells :
#                 'DendriticAlpha' = D1_alpha, 
#                 'NeutrophilAlpha' = D2_alpha,
#                 'allDendriticCommon' = D_alpha)
# GOplatelets <- list(# Platelets : 
#                 'PlateletsCommon' = P_common, 'PlateletsAlpha' = P_alpha, 'PlatetletsLambda' = P_lambda)
GOTcells <- list(# T cells : 
                'Tcommon' = T1_common,'TAlpha' = T1_alpha, 'TLambda' = T1_lambda,
                'CD4hcommon' = T2_common,'CD4hAlpha' = T2_alpha, 'CD4hLambda' = T2_lambda,
                'naiveCD8common' = T3_common,'naiveCD8alpha' = T3_alpha,'naiveCD8Lambda' = T3_lambda,
                'NKTcommon' = T4_common, 'NKTAlpha' = T4_alpha, 'NKTLambda' = T4_lambda, 
                'cCD8common' = T5_common, 'cCD8Alpha' = T5_alpha, 'cCD8Lambda' = T5_lambda,
                'Tregscommon' = T6_common, 'TregsAlpha' = T6_alpha, 'TregsLambda' = T6_lambda,
                # 'CD4Alpha' = T7_alpha, 'CD4Alpha' = T7_alpha,'CD4Alpha' = T7_alp,
               'NKcommon' = T8_common, 'NKAlpha' = T8_alpha,'NKLambda' = T8_lambda,
                'AllTcommon' = T_common,'AllTalpha' = T_alpha,'AllTlambda' = T_lambda)

GOBcells <- list('Bcommon' = B_common, 'BAlpha' = B_alpha, 'BLambda' = B_lambda)
                
openxlsx::write.xlsx(GOmyeloid, file = "honours/results/FinalIndex/GOAnalysis/GOmyeloid.xlsx") 
#openxlsx::write.xlsx(GOdendritic, file = "honours/results/FinalIndex/GOAnalysis/GOdendritic.xlsx") 
#openxlsx::write.xlsx(GOplatelets, file = "honours/results/FinalIndex/GOAnalysis/GOplatelets.xlsx")
openxlsx::write.xlsx(GOTcells, file = "honours/results/FinalIndex/GOAnalysis/GOTcells.xlsx")
openxlsx::write.xlsx(GOBcells, file = "honours/results/FinalIndex/GOAnalysis/GOBcells.xlsx")

##### [4] Gage KEGG analysis #####

library(pathview)
library(gage)
library(gageData)

# reference data sets from gageData : 

data("go.sets.hs") # GO analysis
data("go.subs.hs") # GO analysis 
gobpsets = go.sets.hs[go.subs.hs$BP] # GO analysis 

data(kegg.sets.hs) # KEGG analysis 
data(sigmet.idx.hs) # subset of KEGG sets with only signalling & metabolic pathways 
# kegg.sets.hs = kegg.sets.hs[sigment.idx.hs] # subsetting KEGG sets to only include metabolic processes 


# creating data with log2fold changes and entrex gene IDs :

M1Afoldchanges = M1AlphaResponse$avg_log2FC
names(M1Afoldchanges) = M1AlphaResponse$entrezgene_id
M1Lfoldchanges = M1LambdaResponse$avg_log2FC
names(M1Lfoldchanges) = M1LambdaResponse$entrezgene_id

M2Afoldchanges = M2AlphaResponse$avg_log2FC
names(M2Afoldchanges) = M2AlphaResponse$entrezgene_id
M2Lfoldchanges = M2LambdaResponse$avg_log2FC
names(M2Lfoldchanges) = M2LambdaResponse$entrezgene_id

L1Afoldchanges = L1AlphaResponse$avg_log2FC
names(L1Afoldchanges) = L1AlphaResponse$entrezgene_id
L1Lfoldchanges = L1LambdaResponse$avg_log2FC
names(L1Lfoldchanges) = L1LambdaResponse$entrezgene_id

L2Afoldchanges = L2AlphaResponse$avg_log2FC
names(L2Afoldchanges) = L2AlphaResponse$entrezgene_id
L2Lfoldchanges = L2LambdaResponse$avg_log2FC
names(L2Lfoldchanges) = L2LambdaResponse$entrezgene_id

L3Afoldchanges = L3AlphaResponse$avg_log2FC
names(L3Afoldchanges) = L3AlphaResponse$entrezgene_id
L3Lfoldchanges = L3LambdaResponse$avg_log2FC
names(L3Lfoldchanges) = L3LambdaResponse$entrezgene_id

L4Afoldchanges = L4AlphaResponse$avg_log2FC
names(L4Afoldchanges) = L4AlphaResponse$entrezgene_id
L4Lfoldchanges = L4LambdaResponse$avg_log2FC
names(L4Lfoldchanges) = L4LambdaResponse$entrezgene_id

L5Afoldchanges = L5AlphaResponse$avg_log2FC
names(L5Afoldchanges) = L5AlphaResponse$entrezgene_id
L5Lfoldchanges = L5LambdaResponse$avg_log2FC
names(L5Lfoldchanges) = L5LambdaResponse$entrezgene_id

L6Afoldchanges = L6AlphaResponse$avg_log2FC
names(L6Afoldchanges) = L6AlphaResponse$entrezgene_id
L6Lfoldchanges = L6LambdaResponse$avg_log2FC
names(L6Lfoldchanges) = L6LambdaResponse$entrezgene_id


# GO Biological processes for DEGs using gage :(already done this with clusterprofiler)

M1Agobp = gage(exprs = M1Afoldchanges, # DEGs from our dataset
               gsets = gobpsets, # gageData sets subsetted for BP 
               same.dir = TRUE) # looks in same direction for up & down regulated genes)

# KEGG pathways for DEGs using gage :

# M1 Alpha
M1Akegg = gage(exprs = M1Afoldchanges, # DEGs from our dataset
               gsets = kegg.sets.hs, # gageData sets
               same.dir = TRUE)
saveRDS(M1Akegg, "M1Alpha/M1Akegg.rds")
M1Akeggpathways = data.frame(id = rownames(M1Akegg$greater), M1Akegg$greater) %>% #list with rownames and subset out for name in front
                               tibble::as_tibble() %>% 
                              filter(row_number() <= 20) %>% # top 20 
                                .$id %>%
                                  as.character()
saveRDS(M1Akeggpathways, "M1Alpha/M1Akeggpathways.rds")
M1Akeggids = substr(M1Akeggpathways, start = 1, stop = 8) # just first view characters 
setwd("C:/Users/alice/honours/results/KEGG/") # Make the pathway plots : 
tmp = sapply(M1Akeggids, 
             function(pid) pathview(gene.data = M1Afoldchanges, pathway.id = pid, species = "hsa"))
pathview(gene.data = M1Afoldchanges, pathway.id = "hsa04630", species = "hsa") 


# M1 Lambda 
M1Lkegg = gage(exprs = M1Lfoldchanges, # DEGs from our dataset
               gsets = kegg.sets.hs, # gageData sets
               same.dir = TRUE)
saveRDS(M1Akegg, "M1Lambda/M1Lkegg.rds")
M1Lkeggpathways = data.frame(id = rownames(M1Lkegg$greater), M1Lkegg$greater) %>% #list with rownames and subset out for name in front
  tibble::as_tibble() %>% 
  filter(row_number() <= 20) %>% # top 20 
  .$id %>%
  as.character()
saveRDS(M1Akeggpathways, "M1Lambda/M1Lkeggpathways.rds")
M1Lkeggids = substr(M1Lkeggpathways, start = 1, stop = 8)
tmp = sapply(M1Lkeggids, 
             function(pid) pathview(gene.data = M1Lfoldchanges, pathway.id = pid, species = "hsa"))


# M2 Alpha 
M2Akegg = gage(exprs = M2Afoldchanges, # DEGs from our dataset
               gsets = kegg.sets.hs, # gageData sets
               same.dir = TRUE)
saveRDS(M2Akegg, "M2Alpha/M2Akegg.rds")
M2Akeggpathways = data.frame(id = rownames(M2Akegg$greater), M2Akegg$greater) %>% #list with rownames and subset out for name in front
  tibble::as_tibble() %>% 
  filter(row_number() <= 20) %>% # top 20 
  .$id %>%
  as.character()
saveRDS(M2Akeggpathways, "M2Alpha/M2Akeggpathways.rds")
M2Akeggids = substr(M2Akeggpathways, start = 1, stop = 8) # just first view characters 
setwd("C:/Users/alice/honours/results/KEGG/M2Alpha/") # Make the pathway plots : 
tmp = sapply(M2Akeggids, 
             function(pid) pathview(gene.data = M2Afoldchanges, pathway.id = pid, species = "hsa"))

# M2 Lambda 
M2Lkegg = gage(exprs = M2Lfoldchanges, # DEGs from our dataset
               gsets = kegg.sets.hs, # gageData sets
               same.dir = TRUE)
saveRDS(M2Lkegg, "M2Lambda/M2Lkegg.rds")
M2Lkeggpathways = data.frame(id = rownames(M2Lkegg$greater), M2Lkegg$greater) %>% #list with rownames and subset out for name in front
  tibble::as_tibble() %>% 
  filter(row_number() <= 20) %>% # top 20 
  .$id %>%
  as.character()
saveRDS(M2Lkeggpathways, "M2Lambda/M2Lkeggpathways.rds")
M2Lkeggids = substr(M2Lkeggpathways, start = 1, stop = 8)
tmp = sapply(M2Lkeggids, 
             function(pid) pathview(gene.data = M2Lfoldchanges, pathway.id = pid, species = "hsa"))


# Make the pathway plots of specific pathways: 
setwd("C:/Users/alice/honours/results/KEGG/")
pathview(gene.data = M1Afoldchanges, pathway.id = "hsa04630", species = "hsa") 






# emapplot(gse)
# 4 : view output 
require(DOSE)
barplot(gse, showCategory = 20)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")

emapplot(gse)
library(readr)
full_genes <- read_tsv("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/features_gene_names.tsv")

ego <- enrichGO(gene = gene_list, 
         universe = full_genes,
         OrgDb = organism, 
         keyType = "SYMBOL",
         ont = "CC", 
         pAdjustMethod = "BH",
         pvalueCutoff = 0.01, 
         qvalueCutoff = 0.05, 
         readable = TRUE)

ego2 <- pairwise_termsim(gseGO, method = "Wang", semData )


##### checking DEGs 
getwd()
read_csv("honours/results/FinalIndex/DEAnalysis/BAlphaResponse.csv")
read_csv("honours/results/FinalIndex/DEAnalysis/BLambdaResponse.csv")
read_csv("honours/results/FinalIndex/DEAnalysis/TAlphaResponse.csv")
read_csv("honours/results/FinalIndex/DEAnalysis/TLambdaResponse.csv", spec(), show_col_types = TRUE)



