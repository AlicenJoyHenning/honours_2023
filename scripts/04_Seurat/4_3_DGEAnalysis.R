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
library(biomaRt) # for entrezgenes 
library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(patchwork)
library(openxlsx)
library(mart)

TreatmentAnnotated <- readRDS("honours/results/FinalIndex/TreatmentAnnotated.rds")
Idents(TreatmentAnnotated)

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
# Creating metadata column to hold cell type AND stimulation information (as opposed to just cell type) : 

TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$treatment, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the idents to that column 

# Now we can easily refer to a cell & treatment type, to see the options : 
levels(TreatmentAnnotated) # list the options to perform DE on :

# Use FindMarkers() to find the genes that are different between stimulated and untreated cell types
?FindMarkers()

# MYELOID CELL TYPES : 
# 1. Monocytes 
system.time(M1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "monocytes_alpha", ident.2 = "monocytes_untreated", sep = "_"))
M1AlphaResponse$gene <- rownames(M1AlphaResponse) # puts gene names into a column
M1AlphaResponse <- filter(M1AlphaResponse, p_val_adj < 0.05)
down <- filter(M1AlphaResponse, avg_log2FC < 0)
up <- filter(M1AlphaResponse, avg_log2FC >= 0)
MonoAlphaResponse <- data.frame(Gene = M1AlphaResponse$gene, Log2FoldChange = M1AlphaResponse$avg_log2FC)  # for cluster profiler 
write.csv(MonoAlphaResponse, "honours/results/FinalIndex/DEAnalysis/monocytesAlphaResponse.csv", row.names = FALSE)

system.time(M1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "monocytes_lambda", ident.2 = "monocytes_untreated", sep = "_")) 
M1LambdaResponse$gene <- rownames(M1LambdaResponse) 
M1LambdaResponse <- filter(M1LambdaResponse, p_val_adj < 0.05)
down <- filter(M1LambdaResponse, avg_log2FC < 0)
up <- filter(M1LambdaResponse, avg_log2FC >= 0)
MonoLambdaResponse <- data.frame(Gene = M1LambdaResponse$gene, Log2FoldChange = M1LambdaResponse$avg_log2FC)  
write.csv(MonoLambdaResponse, "honours/results/FinalIndex/DEAnalysis/monocytesLambdaResponse.csv", row.names = FALSE)


# 2. Neutrophils 
system.time(M2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_alpha", ident.2 = "neutrophils_untreated", sep = "_"))  
M2AlphaResponse$gene <- rownames(M2AlphaResponse)
M2AlphaResponse <- filter(M2AlphaResponse, p_val_adj < 0.05)
down <- filter(M2AlphaResponse, avg_log2FC < 0)
up <- filter(M2AlphaResponse, avg_log2FC >= 0)
neuAlphaResponse <- data.frame(Gene = M2AlphaResponse$gene, Log2FoldChange = M2AlphaResponse$avg_log2FC)  
write.csv(neuAlphaResponse, "honours/results/FinalIndex/DEAnalysis/neutrophilsAlphaResponse.csv", row.names = FALSE)
# size = 1507 DEGs 

system.time(M2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_lambda", ident.2 = "neutrophils_untreated", sep = "_")) 
M2LambdaResponse$gene <- rownames(M2LambdaResponse) 
M2LambdaResponse <- filter(M2LambdaResponse, p_val_adj < 0.05)
down <- filter(M2LambdaResponse, avg_log2FC < 0)
up <- filter(M2LambdaResponse, avg_log2FC >= 0)
M2LambdaResponse <- data.frame(Gene = M2LambdaResponse$gene, Log2FoldChange = M2LambdaResponse$avg_log2FC)  
write.csv(M2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/neutrophilsLambdaResponse.csv", row.names = FALSE)

# DENDRITIC CELL TYPES 
# 1. DCs 
system.time(D1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_alpha", ident.2 = "DCs_untreated", sep = "_"))
D1AlphaResponse$gene <- rownames(D1AlphaResponse)
D1AlphaResponse <- filter(D1AlphaResponse, p_val_adj < 0.05)
down <- filter(D1AlphaResponse, avg_log2FC < 0)
up <- filter(D1AlphaResponse, avg_log2FC >= 0)
D1AlphaResponse <- data.frame(Gene = D1AlphaResponse$gene, Log2FoldChange = D1AlphaResponse$avg_log2FC)  
write.csv(D1AlphaResponse, "honours/results/FinalIndex/DEAnalysis/dcsAlphaResponse.csv", row.names = FALSE)

system.time(D1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_lambda", ident.2 = "DCs_untreated", sep = "_")) 
D1LambdaResponse$gene <- rownames(D1LambdaResponse) 
D1LambdaResponse <- filter(D1LambdaResponse, p_val_adj < 0.05)
down <- filter(D1LambdaResponse, avg_log2FC < 0)
up <- filter(D1LambdaResponse, avg_log2FC >= 0)
D1LambdaResponse <- data.frame(Gene = D1LambdaResponse$gene, Log2FoldChange = D1LambdaResponse$avg_log2FC)  
write.csv(D1LambdaResponse, "honours/results/FinalIndex/DEAnalysis/dcsLambdaResponse.csv", row.names = FALSE)

#2. mDCs
system.time(D2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "mDCs_alpha", ident.2 = "mDCs_untreated", sep = "_")) 
D2AlphaResponse$gene <- rownames(D2AlphaResponse)
D2AlphaResponse <- filter(D2AlphaResponse, p_val_adj < 0.05)
down <- filter(D2AlphaResponse, avg_log2FC < 0)
up <- filter(D2AlphaResponse, avg_log2FC >= 0)
D2AlphaResponse <- data.frame(Gene = D2AlphaResponse$gene, Log2FoldChange = D2AlphaResponse$avg_log2FC)  
write.csv(D2AlphaResponse, "honours/results/FinalIndex/DEAnalysis/mdcsAlphaResponse.csv", row.names = FALSE)

system.time(D2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "mDCs_lambda", ident.2 = "mDCs_untreated", sep = "_")) 
D2LambdaResponse$gene <- rownames(D2LambdaResponse) 
D2LambdaResponse <- filter(D2LambdaResponse, p_val_adj < 0.05)
down <- filter(D2LambdaResponse, avg_log2FC < 0)
up <- filter(D2LambdaResponse, avg_log2FC >= 0)
D2LambdaResponse <- data.frame(Gene = D2LambdaResponse$gene, Log2FoldChange = D2LambdaResponse$avg_log2FC)  
write.csv(D2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/mdcsLambdaResponse.csv", row.names = FALSE)

# 3. pDCs : only in alpha dataset, therefore no DEG
# system.time(D3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "pDCs_alpha", ident.2 = "pDCs_untreated", sep = "_")) 
# D3AlphaResponse$gene <- rownames(D3AlphaResponse)
# D3AlphaResponse <- filter(D3AlphaResponse, p_val_adj < 0.05)
# down <- filter(D3AlphaResponse, avg_log2FC < 0)
# up <- filter(D3AlphaResponse, avg_log2FC >= 0)
# D3AlphaResponse <- data.frame(Gene = D3AlphaResponse$gene, Log2FoldChange = D3AlphaResponse$avg_log2FC)  
# write.csv(D3AlphaResponse, "honours/results/FinalIndex/DEAnalysis/pdcsAlphaResponse.csv", row.names = FALSE)
# 
# system.time(D3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "pDCs_lambda", ident.2 = "pDCs_untreated", sep = "_")) 
# D3LambdaResponse$gene <- rownames(D3LambdaResponse) 
# D3LambdaResponse <- filter(D3LambdaResponse, p_val_adj < 0.05)
# down <- filter(D3LambdaResponse, avg_log2FC < 0)
# up <- filter(D3LambdaResponse, avg_log2FC >= 0)
# D3LambdaResponse <- data.frame(Gene = D3LambdaResponse$gene, Log2FoldChange = D3LambdaResponse$avg_log2FC)  
# write.csv(D3LambdaResponse, "honours/results/FinalIndex/DEAnalysis/pdcsLambdaResponse.csv", row.names = FALSE)

# PLATELETS : 
system.time(PAlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_alpha", ident.2 = "platelets_untreated", sep = "_"))
PAlphaResponse$gene <- rownames(PAlphaResponse)
PAlphaResponse <- filter(PAlphaResponse, p_val_adj < 0.05)
down <- filter(PAlphaResponse, avg_log2FC < 0)
up <- filter(PAlphaResponse, avg_log2FC >= 0)
PAlphaResponse <- data.frame(Gene = PAlphaResponse$gene, Log2FoldChange = PAlphaResponse$avg_log2FC)  
write.csv(PAlphaResponse, "honours/results/FinalIndex/DEAnalysis/plateletsAlphaResponse.csv", row.names = FALSE)

system.time(PLambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_lambda", ident.2 = "platelets_untreated", sep = "_"))
PLambdaResponse$gene <- rownames(PLambdaResponse) 
PLambdaResponse <- filter(PLambdaResponse, p_val_adj < 0.05)
down <- filter(PLambdaResponse, avg_log2FC < 0)
up <- filter(PLambdaResponse, avg_log2FC >= 0)
PLambdaResponse <- data.frame(Gene = PLambdaResponse$gene, Log2FoldChange = PLambdaResponse$avg_log2FC)  
write.csv(PLambdaResponse, "honours/results/FinalIndex/DEAnalysis/plateletsLambdaResponse.csv", row.names = FALSE)

# T CELL TYPES : 
# 1. T cells
system.time(T1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_alpha", ident.2 = "T_untreated", sep = "_")) 
T1AlphaResponse$gene <- rownames(T1AlphaResponse)
T1AlphaResponse <- filter(T1AlphaResponse, p_val_adj < 0.05)
down <- filter(T1AlphaResponse, avg_log2FC < 0)
up <- filter(T1AlphaResponse, avg_log2FC >= 0)
T1AlphaResponse <- data.frame(Gene = T1AlphaResponse$gene, Log2FoldChange = T1AlphaResponse$avg_log2FC)  
write.csv(T1AlphaResponse, "honours/results/FinalIndex/DEAnalysis/TAlphaResponse.csv", row.names = FALSE)

system.time(T1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_lambda", ident.2 = "T_untreated", sep = "_")) 
T1LambdaResponse$gene <- rownames(T1LambdaResponse) 
T1LambdaResponse <- filter(T1LambdaResponse, p_val_adj < 0.05)
down <- filter(T1LambdaResponse, avg_log2FC < 0)
up <- filter(T1LambdaResponse, avg_log2FC >= 0)
T1LambdaResponse <- data.frame(Gene = T1LambdaResponse$gene, Log2FoldChange = T1LambdaResponse$avg_log2FC)  
write.csv(T1LambdaResponse, "honours/results/FinalIndex/DEAnalysis/TLambdaResponse.csv", row.names = FALSE)

# 2. CD4+ helper T 
system.time(T2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_helper_alpha", ident.2 = "CD4+_helper_untreated", sep = "_"))  
T2AlphaResponse$gene <- rownames(T2AlphaResponse) 
T2AlphaResponse  <- filter(T2AlphaResponse , p_val_adj < 0.05)
down <- filter(T2AlphaResponse, avg_log2FC < 0)
up <- filter(T2AlphaResponse, avg_log2FC >= 0)
T2AlphaResponse <- data.frame(Gene = T2AlphaResponse$gene, Log2FoldChange = T2AlphaResponse$avg_log2FC)  
write.csv(T2AlphaResponse, "honours/results/FinalIndex/DEAnalysis/CD4hAlphaResponse.csv", row.names = FALSE)

system.time(T2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_helper_lambda", ident.2 = "CD4+_helper_untreated", sep = "_"))  
T2LambdaResponse$gene <- rownames(T2LambdaResponse)
T2LambdaResponse <- filter(T2LambdaResponse, p_val_adj < 0.05)
down <- filter(T2LambdaResponse, avg_log2FC < 0)
up <- filter(T2LambdaResponse, avg_log2FC >= 0)
T2LambdaResponse <- data.frame(Gene = T2LambdaResponse$gene, Log2FoldChange = T2LambdaResponse$avg_log2FC)  
write.csv(T2LambdaResponse, "honours/results/FinalIndex/DEAnalysis/CD4hLambdaResponse.csv", row.names = FALSE)

# 3. Naive CD8 T cells 
system.time(T3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "naive_CD8+_T_alpha", ident.2 = "naive_CD8+_T_untreated", sep = "_")) 
T3AlphaResponse$gene <- rownames(T3AlphaResponse) 
T3AlphaResponse <- filter(T3AlphaResponse, p_val_adj < 0.05)
down <- filter(T3AlphaResponse, avg_log2FC < 0)
up <- filter(T3AlphaResponse, avg_log2FC >= 0)
T3AlphaResponse <- data.frame(Gene = T3AlphaResponse$gene, Log2FoldChange = T3AlphaResponse$avg_log2FC)  
write.csv(T3AlphaResponse, "honours/results/FinalIndex/DEAnalysis/naivecd8AlphaResponse.csv", row.names = FALSE)

system.time(T3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "naive_CD8+_T_lambda", ident.2 = "naive_CD8+_T_untreated", sep = "_"))
T3LambdaResponse$gene <- rownames(T3LambdaResponse)
T3LambdaResponse <- filter(T3LambdaResponse, p_val_adj < 0.05)
down <- filter(T3LambdaResponse, avg_log2FC < 0)
up <- filter(T3LambdaResponse, avg_log2FC >= 0)
T3LambdaResponse <- data.frame(Gene = T3LambdaResponse$gene, Log2FoldChange = T3LambdaResponse$avg_log2FC)  
write.csv(T3LambdaResponse, "honours/results/FinalIndex/DEAnalysis/naivecd8LambdaResponse.csv", row.names = FALSE)
 
# 4. NKT 
system.time(T4AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NKT_alpha", ident.2 = "NKT_untreated", sep = "_"))  
T4AlphaResponse$gene <- rownames(T4AlphaResponse) 
T4AlphaResponse <- filter(T4AlphaResponse, p_val_adj < 0.05)
down <- filter(T4AlphaResponse, avg_log2FC < 0)
up <- filter(T4AlphaResponse, avg_log2FC >= 0)
T4AlphaResponse <- data.frame(Gene = T4AlphaResponse$gene, Log2FoldChange = T4AlphaResponse$avg_log2FC)  
write.csv(T4AlphaResponse, "honours/results/FinalIndex/DEAnalysis/NKTAlphaResponse.csv", row.names = FALSE)

system.time(T4LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NKT_lambda", ident.2 = "NKT_untreated", sep = "_"))
T4LambdaResponse$gene <- rownames(T4LambdaResponse)
T4LambdaResponse <- filter(T4LambdaResponse, p_val_adj < 0.05)
down <- filter(T4LambdaResponse, avg_log2FC < 0)
up <- filter(T4LambdaResponse, avg_log2FC >= 0)
T4LambdaResponse <- data.frame(Gene = T4LambdaResponse$gene, Log2FoldChange = T4LambdaResponse$avg_log2FC)  
write.csv(T4LambdaResponse, "honours/results/FinalIndex/DEAnalysis/NKTLambdaResponse.csv", row.names = FALSE)

# 5. cytotoxic cd8 t cells : 
system.time(T5AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cytotoxic_CD8+_T_alpha", ident.2 = "cytotoxic_CD8+_T_untreated", sep = "_")) 
T5AlphaResponse$gene <- rownames(T5AlphaResponse) 
T5AlphaResponse <- filter(T5AlphaResponse, p_val_adj < 0.05)
down <- filter(T5AlphaResponse, avg_log2FC < 0)
up <- filter(T5AlphaResponse, avg_log2FC >= 0)
T5AlphaResponse <- data.frame(Gene = T5AlphaResponse$gene, Log2FoldChange = T5AlphaResponse$avg_log2FC)  
write.csv(T5AlphaResponse, "honours/results/FinalIndex/DEAnalysis/cytotoxiccd8AlphaResponse.csv", row.names = FALSE)

system.time(T5LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cytotoxic_CD8+_T_lambda", ident.2 = "cytotoxic_CD8+_T_untreated", sep = "_")) 
T5LambdaResponse$gene <- rownames(T5LambdaResponse)
T5LambdaResponse <- filter(T5LambdaResponse, p_val_adj < 0.05)
down <- filter(T5LambdaResponse, avg_log2FC < 0)
up <- filter(T5LambdaResponse, avg_log2FC >= 0)
T5LambdaResponse <- data.frame(Gene = T5LambdaResponse$gene, Log2FoldChange = T5LambdaResponse$avg_log2FC)  
write.csv(T5LambdaResponse, "honours/results/FinalIndex/DEAnalysis/cytotoxiccd8LambdaResponse.csv", row.names = FALSE)

# 6. Tregs : 
system.time(T6AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_alpha", ident.2 = "Tregs_untreated", sep = "_")) 
T6AlphaResponse$gene <- rownames(T6AlphaResponse) 
T6AlphaResponse <- filter(T6AlphaResponse, p_val_adj < 0.05)
down <- filter(T6AlphaResponse, avg_log2FC < 0)
up <- filter(T6AlphaResponse, avg_log2FC >= 0)
T6AlphaResponse <- data.frame(Gene = T6AlphaResponse$gene, Log2FoldChange = T6AlphaResponse$avg_log2FC)  
write.csv(T6AlphaResponse, "honours/results/FinalIndex/DEAnalysis/tregsAlphaResponse.csv", row.names = FALSE)

system.time(T6LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_lambda", ident.2 = "Tregs_untreated", sep = "_")) 
T6LambdaResponse$gene <- rownames(T6LambdaResponse)
T6LambdaResponse <- filter(T6LambdaResponse, p_val_adj < 0.05)
down <- filter(T6LambdaResponse, avg_log2FC < 0)
up <- filter(T6LambdaResponse, avg_log2FC >= 0)
T6LambdaResponse <- data.frame(Gene = T6LambdaResponse$gene, Log2FoldChange = T6LambdaResponse$avg_log2FC)  
write.csv(T6LambdaResponse, "honours/results/FinalIndex/DEAnalysis/tregsLambdaResponse.csv", row.names = FALSE)
 
# 7. CD4+ T cells : 
system.time(T7AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_T_alpha", ident.2 = "CD4+_T_untreated", sep = "_")) 
T7AlphaResponse$gene <- rownames(T7AlphaResponse) 
T7AlphaResponse <- filter(T7AlphaResponse, p_val_adj < 0.05)
down <- filter(T7AlphaResponse, avg_log2FC < 0)
up <- filter(T7AlphaResponse, avg_log2FC >= 0)
T7AlphaResponse <- data.frame(Gene = T7AlphaResponse$gene, Log2FoldChange = T7AlphaResponse$avg_log2FC)  
write.csv(T7AlphaResponse, "honours/results/FinalIndex/DEAnalysis/cd4AlphaResponse.csv", row.names = FALSE)

system.time(T7LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4+_T_lambda", ident.2 = "CD4+_T_untreated", sep = "_"))
T7LambdaResponse$gene <- rownames(T7LambdaResponse)
T7LambdaResponse <- filter(T7LambdaResponse, p_val_adj < 0.05)
down <- filter(T7LambdaResponse, avg_log2FC < 0)
up <- filter(T7LambdaResponse, avg_log2FC >= 0)
T7LambdaResponse <- data.frame(Gene = T7LambdaResponse$gene, Log2FoldChange = T7LambdaResponse$avg_log2FC)  
write.csv(T7LambdaResponse, "honours/results/FinalIndex/DEAnalysis/cd4LambdaResponse.csv", row.names = FALSE)

# 8. NK cells 
system.time(T8AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_alpha", ident.2 = "NK_untreated", sep = "_"))
T8AlphaResponse$gene <- rownames(T8AlphaResponse) 
T8AlphaResponse <- filter(T8AlphaResponse, p_val_adj < 0.05)
down <- filter(T8AlphaResponse, avg_log2FC < 0)
up <- filter(T8AlphaResponse, avg_log2FC >= 0)
T8AlphaResponse <- data.frame(Gene = T8AlphaResponse$gene, Log2FoldChange = T8AlphaResponse$avg_log2FC)  
write.csv(T8AlphaResponse, "honours/results/FinalIndex/DEAnalysis/NKAlphaResponse.csv", row.names = FALSE)

system.time(T8LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_lambda", ident.2 = "NK_untreated", sep = "_")) 
T8LambdaResponse$gene <- rownames(T8LambdaResponse)
T8LambdaResponse <- filter(T8LambdaResponse, p_val_adj < 0.05)
down <- filter(T8LambdaResponse, avg_log2FC < 0)
up <- filter(T8LambdaResponse, avg_log2FC >= 0)
T8LambdaResponse <- data.frame(Gene = T8LambdaResponse$gene, Log2FoldChange = T8LambdaResponse$avg_log2FC)  
write.csv(T8LambdaResponse, "honours/results/FinalIndex/DEAnalysis/NKLambdaResponse.csv", row.names = FALSE)

# B CELLS 
system.time(BAlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_alpha", ident.2 = "B_untreated", sep = "_")) 
BAlphaResponse$gene <- rownames(BAlphaResponse) 
BAlphaResponse <- filter(BAlphaResponse, p_val_adj < 0.05)
down <- filter(BAlphaResponse, avg_log2FC < 0)
up <- filter(BAlphaResponse, avg_log2FC >= 0)
BAlphaResponse <- data.frame(Gene = BAlphaResponse$gene, Log2FoldChange = BAlphaResponse$avg_log2FC)  
write.csv(BAlphaResponse, "honours/results/FinalIndex/DEAnalysis/BAlphaResponse.csv", row.names = FALSE)

system.time(BLambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_lambda", ident.2 = "B_untreated", sep = "_"))
BLambdaResponse$gene <- rownames(BLambdaResponse)
BLambdaResponse <- filter(BLambdaResponse, p_val_adj < 0.05)
down <- filter(BLambdaResponse, avg_log2FC < 0)
up <- filter(BLambdaResponse, avg_log2FC >= 0)
BLambdaResponse <- data.frame(Gene = BLambdaResponse$gene, Log2FoldChange = BLambdaResponse$avg_log2FC)  
write.csv(BLambdaResponse, "honours/results/FinalIndex/DEAnalysis/BLambdaResponse.csv", row.names = FALSE)



##### [2.2] DE with FindMarkers() on broad cell types () #####

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
MAResponse <- data.frame(Gene = MAResponse$gene, Log2FoldChange = MAResponse$avg_log2FC)  # for cluster profiler 
write.csv(MAResponse, "honours/results/FinalIndex/DEAnalysis/MyeloidAlphaResponse.csv", row.names = FALSE)

system.time(MLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "myeloid_lambda", ident.2 = "myeloid_untreated", sep = "_"))
MLResponse$gene <- rownames(MLResponse) # puts gene names into a column
MLResponse <- filter(MLResponse, p_val_adj < 0.05)
down <- filter(MLResponse, avg_log2FC < 0)
up <- filter(MLResponse, avg_log2FC >= 0)
MLResponse <- data.frame(Gene = MLResponse$gene, Log2FoldChange = MLResponse$avg_log2FC)  # for cluster profiler 
write.csv(MLResponse, "honours/results/FinalIndex/DEAnalysis/MyeloidLambdaResponse.csv", row.names = FALSE)

# Dendritic cell types : 
system.time(DAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_alpha", ident.2 = "DCs_untreated", sep = "_"))
DAResponse$gene <- rownames(DAResponse) # puts gene names into a column
DAResponse <- filter(DAResponse, p_val_adj < 0.05)
down <- filter(DAResponse, avg_log2FC < 0)
up <- filter(DAResponse, avg_log2FC >= 0)
DAResponse <- data.frame(Gene = DAResponse$gene, Log2FoldChange = DAResponse$avg_log2FC)  # for cluster profiler 
write.csv(DAResponse, "honours/results/FinalIndex/DEAnalysis/DendriticAlphaResponse.csv", row.names = FALSE)

system.time(DLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "DCs_lambda", ident.2 = "DCs_untreated", sep = "_"))
DLResponse$gene <- rownames(DLResponse) # puts gene names into a column
DLResponse <- filter(DLResponse, p_val_adj < 0.05)
down <- filter(DLResponse, avg_log2FC < 0)
up <- filter(DLResponse, avg_log2FC >= 0)
DLResponse <- data.frame(Gene = DLResponse$gene, Log2FoldChange = DLResponse$avg_log2FC)  # for cluster profiler 
write.csv(DLResponse, "honours/results/FinalIndex/DEAnalysis/DendriticLambdaResponse.csv", row.names = FALSE)

# Platelets : 
system.time(PAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_alpha", ident.2 = "platelets_untreated", sep = "_"))
PAResponse$gene <- rownames(PAResponse) # puts gene names into a column
PAResponse <- filter(PAResponse, p_val_adj < 0.05)
down <- filter(PAResponse, avg_log2FC < 0)
up <- filter(PAResponse, avg_log2FC >= 0)
PAResponse <- data.frame(Gene = PAResponse$gene, Log2FoldChange = PAResponse$avg_log2FC)  # for cluster profiler 
write.csv(PAResponse, "honours/results/FinalIndex/DEAnalysis/PlateletsAlphaResponse.csv", row.names = FALSE)

system.time(PLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "platelets_lambda", ident.2 = "platelets_untreated", sep = "_"))
PLResponse$gene <- rownames(PLResponse) # puts gene names into a column
PLResponse <- filter(PLResponse, p_val_adj < 0.05)
down <- filter(PLResponse, avg_log2FC < 0)
up <- filter(PLResponse, avg_log2FC >= 0)
PLResponse <- data.frame(Gene = PLResponse$gene, Log2FoldChange = PLResponse$avg_log2FC)  # for cluster profiler 
write.csv(PLResponse, "honours/results/FinalIndex/DEAnalysis/PlateletsLambdaResponse.csv", row.names = FALSE)

# T cells : 
system.time(TAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_alpha", ident.2 = "T_untreated", sep = "_"))
TAResponse$gene <- rownames(TAResponse) # puts gene names into a column
TAResponse <- filter(TAResponse, p_val_adj < 0.05)
down <- filter(TAResponse, avg_log2FC < 0)
up <- filter(TAResponse, avg_log2FC >= 0)
TAResponse <- data.frame(Gene = TAResponse$gene, Log2FoldChange = TAResponse$avg_log2FC)  # for cluster profiler 
write.csv(TAResponse, "honours/results/FinalIndex/DEAnalysis/TAlphaResponse.csv", row.names = FALSE)

system.time(TLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "T_lambda", ident.2 = "T_untreated", sep = "_"))
TLResponse$gene <- rownames(TLResponse) # puts gene names into a column
TLResponse <- filter(TLResponse, p_val_adj < 0.05)
down <- filter(TLResponse, avg_log2FC < 0)
up <- filter(TLResponse, avg_log2FC >= 0)
TLResponse <- data.frame(Gene = TLResponse$gene, Log2FoldChange = TLResponse$avg_log2FC)  # for cluster profiler 
write.csv(TLResponse, "honours/results/FinalIndex/DEAnalysis/TLambdaResponse.csv", row.names = FALSE)

# B cells : 
system.time(BAResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_alpha", ident.2 = "B_untreated", sep = "_"))
BAResponse$gene <- rownames(BAResponse) # puts gene names into a column
BAResponse <- filter(BAResponse, p_val_adj < 0.05)
down <- filter(BAResponse, avg_log2FC < 0)
up <- filter(BAResponse, avg_log2FC >= 0)
BAResponse <- data.frame(Gene = BAResponse$gene, Log2FoldChange = BAResponse$avg_log2FC)  # for cluster profiler 
write.csv(BAResponse, "honours/results/FinalIndex/DEAnalysis/BAlphaResponse.csv", row.names = FALSE)

system.time(BLResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_lambda", ident.2 = "B_untreated", sep = "_"))
BLResponse$gene <- rownames(BLResponse) # puts gene names into a column
BLResponse <- filter(BLResponse, p_val_adj < 0.05)
down <- filter(BLResponse, avg_log2FC < 0)
up <- filter(BLResponse, avg_log2FC >= 0)
BLResponse <- data.frame(Gene = BLResponse$gene, Log2FoldChange = BLResponse$avg_log2FC)  # for cluster profiler 
write.csv(BLResponse, "honours/results/FinalIndex/DEAnalysis/BLambdaResponse.csv", row.names = FALSE)


##### [2.3] Looking for common & unique DEGs #####

# [a] To find common DEG between 2 treatments compared to the control (same pathways?) : 
B_cells_common_genes <- intersect(BAResponse$Gene, BLResponse$Gene)                # identify common genes
B_cells_common_alpha <- BAResponse[BAResponse$Gene %in% B_cells_common_genes, ]     # Extract common genes from each dataframe 
B_cells_common_lambda <- BLResponse[BLResponse$Gene %in% B_cells_common_genes, ]
B_cells_common_dataframe <- merge(B_cells_common_alpha, B_cells_common_lambda, by = 'Gene', all = TRUE) # # merge together the 2 subseted dataframes to create 1 containing only common genes 

# [b] To find treatment specific DEGs : 
B_cells_alpha_genes <- BAResponse$Gene[!(BAResponse$Gene %in% BLResponse$Gene)] # - (all alpha in lambda) = all alpha not in lambda
B_cells_alpha <- BAResponse[BAResponse$Gene %in% B_cells_alpha_genes, ]


B_cells_lambda_genes <- BLResponse$Gene[!(BLResponse$Gene %in% BAResponse$Gene)]
B_cells_lambda <- BLResponse[BLResponse$Gene %in% B_cells_lambda_genes, ] 

down <- filter(B_cells_alpha, Log2FoldChange < 0)
up <- filter(B_cells_alpha, Log2FoldChange >= 0)


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
listDatasets(mart) # List available datasets for the Mart (214!)

# 3:  Retrieve gene information for human genes :
gene_ID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
                 mart = mart) 

# 4: create a background universe from the 2000 most variable genes 
# Assuming 'mySeurat' is your integrated Seurat object
# Access the integrated assay slot
integrated_assay <- TreatmentAnnotated@assays[["integrated"]]
rna_assay <- rownames(TreatmentAnnotated@assays[['RNA']])
# Get the 2000 most variable genes
top_variable_genes <- rownames(integrated_assay)[order(integrated_assay, decreasing = TRUE)[1:2000]]
universe <- select(org.Hs.eg.db, keys = rna_assay, keytype = "SYMBOL", columns = "ENTREZID")
universe <- universe$ENTREZID




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
BLambdaResponse <- BLambdaResponse[, c(-7, -8)]

# 5 : create object to perform GO analysis on, just the entrezgene  id column : 
# myeloid 
M1AlphaDE <- M1AlphaResponse$entrezgene_id
M1LambdaDE <- M1LambdaResponse$entrezgene_id.x
M2AlphaDE <- M2AlphaResponse$entrezgene_id
M2LambdaDE <- M2LambdaResponse$entrezgene_id
MAlphaDE <- MAResponse$entrezgene_id
MLambdaDE <- MLResponse$entrezgene_id

# dendritic 
D1AlphaDE <- D1AlphaResponse$entrezgene_id
D1LambdaDE <- D1LambdaResponse$entrezgene_id
D2AlphaDE <- D2AlphaResponse$entrezgene_id
D2LambdaDE <- D2LambdaResponse$entrezgene_id
# D3AlphaDE <- D3AlphaResponse$entrezgene_id
# D3LambdaDE <- D3LambdaResponse$entrezgene_id
DAlphaDE <- DAResponse$entrezgene_id
DLambdaDE <- DLResponse$entrezgene_id

# platelets 
PAlphaDE <- PAlphaResponse$entrezgene_id
PLambdaDE <- PLambdaResponse$entrezgene_id

# T cells 
T1AlphaDE <- T1AlphaResponse$entrezgene_id
T1LambdaDE <- T1LambdaResponse$entrezgene_id
T2AlphaDE <- T2AlphaResponse$entrezgene_id
T2LambdaDE <- T2LambdaResponse$entrezgene_id
T3AlphaDE <- T3AlphaResponse$entrezgene_id
T3LambdaDE <- T3LambdaResponse$entrezgene_id
T4AlphaDE <- T4AlphaResponse$entrezgene_id
T4LambdaDE <- T4LambdaResponse$entrezgene_id
T5AlphaDE <- T5AlphaResponse$entrezgene_id
T5LambdaDE <- T5LambdaResponse$entrezgene_id
T6AlphaDE <- T6AlphaResponse$entrezgene_id
T6LambdaDE <- T6LambdaResponse$entrezgene_id
T7AlphaDE <- T7AlphaResponse$entrezgene_id
T7LambdaDE <- T7LambdaResponse$entrezgene_id
T8AlphaDE <- T8AlphaResponse$entrezgene_id
T8LambdaDE <- T8LambdaResponse$entrezgene_id
TAlphaDE <- TAResponse$entrezgene_id
TLambdaDE <- TLResponse$entrezgene_id

# B cells 
BAlphaDE <- BAlphaResponse$entrezgene_id
BLambdaDE <- BLambdaResponse$entrezgene_id

# 6 : perform GO analysis 
?enrichGO


# Myeloid  
M1AlphaEGO <- enrichGO(gene = M1AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH", universe = universe, pvalueCutoff = 0.05)
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
D1LambdaEGO <- enrichGO(gene = D1LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
D1LambdaEGO <- filter(D1LambdaEGO, D1LambdaEGO@result$p.adjust < 0.05)
D2AlphaEGO <- enrichGO(gene = D2AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
D2AlphaEGO <- filter(D2AlphaEGO, D2AlphaEGO@result$p.adjust < 0.05)
D2LambdaEGO <- enrichGO(gene = D2LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
D2LambdaEGO <- filter(D2LambdaEGO, D2LambdaEGO@result$p.adjust < 0.05)
# D3AlphaEGO <- enrichGO(gene = D3AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D3AlphaEGO <- filter(D3AlphaEGO, D3AlphaEGO@result$p.adjust < 0.05)
# D3LambdaEGO <- enrichGO(gene = D3LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# D3LambdaEGO <- filter(D3LambdaEGO, D3LambdaEGO@result$p.adjust < 0.05)
DAlphaEGO <- enrichGO(gene = DAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
DAlphaEGO <- filter(DAlphaEGO, DAlphaEGO@result$p.adjust < 0.05)
DLambdaEGO <- enrichGO(gene = DLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
DLambdaEGO <- filter(DLambdaEGO, DLambdaEGO@result$p.adjust < 0.05)

# Platelets 
PAlphaEGO <- enrichGO(gene = PAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
PAlphaEGO <- filter(PAlphaEGO, PAlphaEGO@result$p.adjust < 0.05)
PLambdaEGO <- enrichGO(gene = PLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
PLambdaEGO <- filter(PLambdaEGO, PLambdaEGO@result$p.adjust < 0.05)

# T cells  
T1AlphaEGO <- enrichGO(gene = T1AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T1AlphaEGO <- filter(T1AlphaEGO, T1AlphaEGO@result$p.adjust < 0.05)
T1_Alpha <- T1AlphaEGO@result$Description
# T1LambdaEGO <- enrichGO(gene = T1LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# T1LambdaEGO <- filter(T1LambdaEGO, T1LambdaEGO@result$p.adjust < 0.05)
T2AlphaEGO <- enrichGO(gene = T2AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T2AlphaEGO <- filter(T2AlphaEGO, T2AlphaEGO@result$p.adjust < 0.05)
T2_Alpha <- T2AlphaEGO@result$Description
# T2LambdaEGO <- enrichGO(gene = T2LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
# T2LambdaEGO <- filter(T2LambdaEGO, T2LambdaEGO@result$p.adjust < 0.05)
T3AlphaEGO <- enrichGO(gene = T3AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T3AlphaEGO <- filter(T3AlphaEGO, T3AlphaEGO@result$p.adjust < 0.05)
T3LambdaEGO <- enrichGO(gene = T3LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T3LambdaEGO <- filter(T3LambdaEGO, T3LambdaEGO@result$p.adjust < 0.05)
T3_Lambda <- T3LambdaEGO@result$Description
T4AlphaEGO <- enrichGO(gene = T4AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T4AlphaEGO <- filter(T4AlphaEGO, T4AlphaEGO@result$p.adjust < 0.05)
T4_Alpha <- T4AlphaEGO@result$Description
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
T7LambdaEGO <- enrichGO(gene = T7LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T7LambdaEGO <- filter(T7LambdaEGO, T7LambdaEGO@result$p.adjust < 0.05)
T8AlphaEGO <- enrichGO(gene = T8AlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T8AlphaEGO <- filter(T8AlphaEGO, T8AlphaEGO@result$p.adjust < 0.05)
T8LambdaEGO <- enrichGO(gene = T8LambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
T8LambdaEGO <- filter(T8LambdaEGO, T8LambdaEGO@result$p.adjust < 0.05)
TAlphaEGO <- enrichGO(gene = TAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
TAlphaEGO <- filter(TAlphaEGO, TAlphaEGO@result$p.adjust < 0.05)
TLambdaEGO <- enrichGO(gene = TLambdaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
TLambdaEGO <- filter(TLambdaEGO, TLambdaEGO@result$p.adjust < 0.05)
T_Alpha <- TAlphaEGO@result$Description

# B CELLS 
BAlphaEGO <- enrichGO(gene = BAlphaDE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.05, universe = universe)
BAlphaEGO <- filter(BAlphaEGO, BAlphaEGO@result$p.adjust < 0.05)
BLambdaEGO <- enrichGO(gene = BLambdaDE, OrgDb = org.Hs.eg.db,keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH",  pvalueCutoff = 0.05, universe = universe)
BLambdaEGO <- filter(BLambdaEGO, BLambdaEGO@result$qvalue < 0.10)


# Then, proceed with enrichGO
clusterProfiler_bp <- enrichGO(BLambdaDE, ont = "BP", OrgDb = org.Hs.eg.db)

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
D1_common <- intersect(D1AlphaEGO@result$Description, D1LambdaEGO@result$Description )                # identify common genes
D1_alpha <- D1AlphaEGO@result$Description[!(D1AlphaEGO@result$Description %in% D1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
D1_lambda <- D1LambdaEGO@result$Description[!(D1LambdaEGO@result$Description %in% D1AlphaEGO@result$Description)]
D1_lambda <- D1LambdaEGO[D1LambdaEGO@result$Description %in% D1_lambda, ] 
D1_lambda_list <- D1_lambda[, c("ID", "p.adjust")] 
write_tsv(D1_lambda_list, "honours/results/FinalIndex/GOAnalysis/D1_lambda_list.tsv")

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

D_common <- intersect(DAlphaEGO@result$Description, DLambdaEGO@result$Description )                # identify common genes
D_alpha <- DAlphaEGO@result$Description[!(DAlphaEGO@result$Description %in% DLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
D_lambda <- DLambdaEGO@result$Description[!(DLambdaEGO@result$Description %in% DAlphaEGO@result$Description)]
D_lambda <- DLambdaEGO[DLambdaEGO@result$Description %in% D_lambda, ] 
D_lambda_list <- D_lambda[, c("ID", "p.adjust")] 
write_tsv(D_lambda_list, "honours/results/FinalIndex/GOAnalysis/D_lambda_list.tsv")

# Platelets : 
P_common <- intersect(PAlphaEGO@result$Description, PLambdaEGO@result$Description )                # identify common genes
P_alpha <- PAlphaEGO@result$Description[!(PAlphaEGO@result$Description %in% PLambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
P_lambda <- PLambdaEGO@result$Description[!(PLambdaEGO@result$Description %in% PAlphaEGO@result$Description)]
P_lambda <- PLambdaEGO[PLambdaEGO@result$Description %in% P_lambda, ] 
P_lambda_list <- P_lambda[, c("ID", "p.adjust")] 
write_tsv(P_lambda_list, "honours/results/FinalIndex/GOAnalysis/P_lambda_list.tsv")

# T cells : 
# T1_common <- intersect(T1AlphaEGO@result$Description, T1LambdaEGO@result$Description )                # identify common genes
T1_alpha <- T1AlphaEGO@result$Description[!(T1AlphaEGO@result$Description %in% T1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# T1_lambda <- T1LambdaEGO@result$Description[!(T1LambdaEGO@result$Description %in% T1AlphaEGO@result$Description)]
# T1_lambda <- T1LambdaEGO[T1LambdaEGO@result$Description %in% T1_lambda, ] 
# T1_lambda_list <- T1_lambda[, c("ID", "p.adjust")] 
# write_tsv(T1_lambda_list, "honours/results/FinalIndex/GOAnalysis/T1_lambda_list.tsv")
# 
# T2_common <- intersect(T2AlphaEGO@result$Description, T2LambdaEGO@result$Description )                # identify common genes
# T2_alpha <- T2AlphaEGO@result$Description[!(T2AlphaEGO@result$Description %in% T2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# T2_lambda <- T2LambdaEGO@result$Description[!(T2LambdaEGO@result$Description %in% T2AlphaEGO@result$Description)]
# T2_lambda <- T2LambdaEGO[T2LambdaEGO@result$Description %in% T2_lambda, ] 
# T2_lambda_list <- T2_lambda[, c("ID", "p.adjust")] 
# write_tsv(T2_lambda_list, "honours/results/FinalIndex/GOAnalysis/T2_lambda_list.tsv")

# T3_common <- intersect(T3AlphaEGO@result$Description, T3LambdaEGO@result$Description )                # identify common genes
# T3_alpha <- T3AlphaEGO@result$Description[!(T3AlphaEGO@result$Description %in% T3LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
# T3_lambda <- T3LambdaEGO@result$Description[!(T3LambdaEGO@result$Description %in% T3AlphaEGO@result$Description)]
# T3_lambda <- T3LambdaEGO[T3LambdaEGO@result$Description %in% T3_lambda, ] 
# T3_lambda_list <- T3_lambda[, c("ID", "p.adjust")] 
# write_tsv(T3_lambda_list, "honours/results/FinalIndex/GOAnalysis/T3_lambda_list.tsv")

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

T7_common <- intersect(T7AlphaEGO@result$Description, T7LambdaEGO@result$Description )                # identify common genes
T7_alpha <- T7AlphaEGO@result$Description[!(T7AlphaEGO@result$Description %in% T7LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T7_lambda <- T7LambdaEGO@result$Description[!(T7LambdaEGO@result$Description %in% T7AlphaEGO@result$Description)]
T7_lambda <- T7LambdaEGO[T7LambdaEGO@result$Description %in% T7_lambda, ] 
T7_lambda_list <- T7_lambda[, c("ID", "p.adjust")] 
write_tsv(T7_lambda_list, "honours/results/FinalIndex/GOAnalysis/T7_lambda_list.tsv")

T8_common <- intersect(T8AlphaEGO@result$Description, T8LambdaEGO@result$Description )                # identify common genes
T8_alpha <- T8AlphaEGO@result$Description[!(T8AlphaEGO@result$Description %in% T8LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
T8_lambda <- T8LambdaEGO@result$Description[!(T8LambdaEGO@result$Description %in% T8AlphaEGO@result$Description)]
T8_lambda <- T8LambdaEGO[T8LambdaEGO@result$Description %in% T8_lambda, ] 
T8_lambda_list <- T8_lambda[, c("ID", "p.adjust")] 
write_tsv(T8_lambda_list, "honours/results/FinalIndex/GOAnalysis/T8_lambda_list.tsv")

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
GOdendritic <- list(# Dendritic cells :
                'DendriticCommon' = D1_common, 'DendriticAlpha' = D1_alpha, 'DendriticLambda' = D1_lambda, 
                'allDendriticCommon' = D_alpha)
GOplatelets <- list(# Platelets : 
                'PlateletsCommon' = P_common, 'PlateletsAlpha' = P_alpha, 'PlatetletsLambda' = P_lambda)
GOTcells <- list(# T cells : 
                'TAlpha' = T1_Alpha,
                'CD4hAlpha' = T2_Alpha, 
                'naiveCD8Lambda' = T3_Lambda,
                'NKTAlpha' = T4_Alpha, 
                'cCD8common' = T5_common, 'cCD8Alpha' = T5_alpha, 'cCD8Lambda' = T5_lambda,
                'Tregscommon' = T6_common, 'TregsAlpha' = T6_alpha, 'TregsLambda' = T6_lambda,
                'CD4Alpha' = T7_alpha, 
               'NKAlpha' = T8_alpha, 
                'AllTalpha' = T_Alpha)
GOBcells <- list('Bcommon' = B_common, 'BAlpha' = B_alpha, 'BLambda' = B_lambda)
                
openxlsx::write.xlsx(GOmyeloid, file = "honours/results/FinalIndex/GOAnalysis/GOmyeloid.xlsx") 
openxlsx::write.xlsx(GOdendritic, file = "honours/results/FinalIndex/GOAnalysis/GOdendritic.xlsx") 
openxlsx::write.xlsx(GOplatelets, file = "honours/results/FinalIndex/GOAnalysis/GOplatelets.xlsx")
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

