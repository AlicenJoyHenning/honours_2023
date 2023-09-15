# DIFFERENTIAL EXPRESSION ANALYSIS (Seurat) 

##### [1.1] Dependencies & load object #####
BiocManager::install("clusterProfiler", version = "3.14")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("gage")
BiocManager::install("gageData")

library(clusterProfiler)
library(pathview)
library(enrichplot)
library(DOSE)
library(biomaRt)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(patchwork)

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



##### [2.2] DE with FIndMarkers() on broad cell types () #####

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


##### [2.2] Extra use of FindMarkers() #####

# [a] To find common DEG between 2 treatments compared to the control (same pathways?) : 
B_cells_common_genes <- intersect(BAResponse$Gene, BLResponse$Gene)                # identify common genes
# [b] To find treatment specific DEGs : 
B_cells_alpha_genes <- BAResponse$gene[!(BAResponse$gene %in% BLResponse$gene)] # - (all alpha in lambda) = all alpha not in lambda
B_cells_alpha <- AlphaResponse[AlphaResponse$gene %in% B_cells_alpha_genes, ]
# size = 1064 genes unique to alpha 

B_cells_lambda_genes <- BLResponse$gene[!(BLResponse$gene %in% BLResponse$gene)]
B_cells_lambda <- LambdaResponse[LambdaResponse$gene %in% B_cells_lambda_genes, ] 
# size = 87 genes unique to lambda
# [c] Store information in excel file 
B_cells <- list('B_cells_Alpha_Response' = AlphaResponse, 'B_cells_Lambda_Response' = LambdaResponse,'B_cells_common_response' = B_cells_common_dataframe, 'B_cells_alpha_specific' = B_cells_alpha,'B_cells_lambda_specific' = B_cells_lambda)
openxlsx::write.xlsx(B_cells, file = "honours/results/DEAnalysis/B_cells.xlsx")

##### [3] ClusterProfiler GO analysis #####

?clusterProfiler # not helpful 
ls("package:clusterProfiler")


# 1 : install and load annotation for desired organism (human) : org.Hs.eg.db
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# 2 : have to convert DE gene list to contain entrezgene_ID BIOMART
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Create a Mart object for the human genome
listDatasets(mart) # List available datasets for the Mart 

# 3:  Retrieve gene information for human genes :
gene_ID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
                 mart = mart)

#### REPITITION STARTS HERE : 
# 4 : merge this with DEGs datasets : 

M1AlphaResponse <- merge(M1AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
M1LambdaResponse <- merge(M1LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

M2AlphaResponse <- merge(M2AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
M2LambdaResponse <- merge(M2LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L1AlphaResponse <- merge(L1AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L1LambdaResponse <- merge(L1LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L2AlphaResponse <- merge(L2AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L2LambdaResponse <- merge(L2LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L3AlphaResponse <- merge(L3AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L3LambdaResponse <- merge(L3LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L4AlphaResponse <- merge(L4AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L4LambdaResponse <- merge(L4LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L5AlphaResponse <- merge(L5AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L5LambdaResponse <- merge(L5LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

L6AlphaResponse <- merge(L6AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")
L6LambdaResponse <- merge(L6LambdaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")

# 5 : create object to perform GO analysis on, just the entrezgene  id column : 

M1AlphaDE <- M1AlphaResponse$entrezgene_id
M1LambdaDE <- M1LambdaResponse$entrezgene_id
M2AlphaDE <- M2AlphaResponse$entrezgene_id
M2LambdaDE <- M2LambdaResponse$entrezgene_id
L1AlphaDE <- L1AlphaResponse$entrezgene_id
L1LambdaDE <- L1LambdaResponse$entrezgene_id
L2AlphaDE <- L2AlphaResponse$entrezgene_id
L2LambdaDE <- L2LambdaResponse$entrezgene_id
L3AlphaDE <- L3AlphaResponse$entrezgene_id
L3LambdaDE <- L3LambdaResponse$entrezgene_id
L4AlphaDE <- L4AlphaResponse$entrezgene_id
L4LambdaDE <- L4LambdaResponse$entrezgene_id
L5AlphaDE <- L5AlphaResponse$entrezgene_id
L5LambdaDE <- L5LambdaResponse$entrezgene_id
L6AlphaDE <- L6AlphaResponse$entrezgene_id
L6LambdaDE <- L6LambdaResponse$entrezgene_id

# 6 : perform GO analysis 
# M1 : Monocytes 
M1AlphaEGO <- enrichGO(gene = M1AlphaDE, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05
                )
M1AlphaEGO <- filter(M1AlphaEGO, M1AlphaEGO@result$p.adjust < 0.05)
M1AlphaList <- M1AlphaEGO@result[, c("ID", "p.adjust")] 
colnames(M1AlphaList) <- c("% GOterm", "enrichment_P-value")
write_tsv(M1AlphaList, "honours/results/DEAnalysis/FortopGO/M1AlphaList.tsv")

M1LambdaEGO <- enrichGO(gene = M1LambdaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
M1LambdaEGO <- filter(M1LambdaEGO, M1LambdaEGO@result$p.adjust < 0.05)

# M2 : Neutrophils

M2AlphaEGO <- enrichGO(gene = M2AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
M2AlphaEGO <- filter(M2AlphaEGO, M2AlphaEGO@result$p.adjust < 0.05)

M2LambdaEGO <- enrichGO(gene = M2LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
M2LambdaEGO <- filter(M2LambdaEGO, M2LambdaEGO@result$p.adjust < 0.05)

# L1 : CD4 helper T cells 
L1AlphaEGO <- enrichGO(gene = L1AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L1AlphaEGO <- filter(L1AlphaEGO, L1AlphaEGO@result$p.adjust < 0.05)

L1LambdaEGO <- enrichGO(gene = L1LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L1LambdaEGO <- filter(L1LambdaEGO, L1LambdaEGO@result$p.adjust < 0.05)


# L2 : CD4 naive T cells 

L2AlphaEGO <- enrichGO(gene = L2AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L2AlphaEGO <- filter(L2AlphaEGO, L2AlphaEGO@result$p.adjust < 0.05)


L2LambdaEGO <- enrichGO(gene = L2LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L2LambdaEGO <- filter(L2LambdaEGO, L2LambdaEGO@result$p.adjust < 0.05)

# L3 : T regulatory cells 

L3AlphaEGO <- enrichGO(gene = L3AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L3AlphaEGO <- filter(L3AlphaEGO, L3AlphaEGO@result$p.adjust < 0.05)

L3LambdaEGO <- enrichGO(gene = L3LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L3LambdaEGO <- filter(L3LambdaEGO, L3LambdaEGO@result$p.adjust < 0.05)

# L4 : CD8 T cells 
L4AlphaEGO <- enrichGO(gene = L4AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L4AlphaEGO <- filter(L4AlphaEGO, L4AlphaEGO@result$p.adjust < 0.05)

L4LambdaEGO <- enrichGO(gene = L4LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L4LambdaEGO <- filter(L4LambdaEGO, L4LambdaEGO@result$p.adjust < 0.05)

# L5 : NK cells 
L5AlphaEGO <- enrichGO(gene = L5AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L5AlphaEGO <- filter(L5AlphaEGO, L5AlphaEGO@result$p.adjust < 0.05)

L5LambdaEGO <- enrichGO(gene = L5LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L5LambdaEGO <- filter(L5LambdaEGO, L5LambdaEGO@result$p.adjust < 0.05)


# L6 : B cells 

L6AlphaEGO <- enrichGO(gene = L6AlphaDE, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05
)
L6AlphaEGO <- filter(L6AlphaEGO, L6AlphaEGO@result$p.adjust < 0.05)

L6LambdaEGO <- enrichGO(gene = L6LambdaDE, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05
)
L6LambdaEGO <- filter(L6LambdaEGO, L6LambdaEGO@result$p.adjust < 0.05)


# 6 : of the descriptions found in alpha & lambda set, subset according to unique to alpha and lambda and the intersection 

L6_common <- intersect(L6AlphaEGO@result$Description, L6LambdaEGO@result$Description )                # identify common genes
L6_alpha <- L6AlphaEGO@result$Description[!(L6AlphaEGO@result$Description %in% L6LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L6_lambda <- L6LambdaEGO@result$Description[!(L6LambdaEGO@result$Description %in% L6AlphaEGO@result$Description)]
L6_lambda <- L6LambdaEGO[L6LambdaEGO@result$Description %in% L6_lambda, ] 
L6_lambda_list <- L6_lambda[, c("ID", "p.adjust")] 
write_tsv(L6_lambda_list, "honours/results/DEAnalysis/FortopGO/L6_lambda_list.tsv")

L5_common <- intersect(L5AlphaEGO@result$Description, L5LambdaEGO@result$Description )                # identify common genes
L5_alpha <- L5AlphaEGO@result$Description[!(L5AlphaEGO@result$Description %in% L5LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L5_lambda <- L5LambdaEGO@result$Description[!(L5LambdaEGO@result$Description %in% L5AlphaEGO@result$Description)]
L5_lambda <- L5LambdaEGO[L5LambdaEGO@result$Description %in% L5_lambda, ] 
L5_lambda_list <- L5_lambda[, c("ID", "p.adjust")] 
write_tsv(L5_lambda_list, "honours/results/DEAnalysis/FortopGO/L5_lambda_list.tsv")

L4_common <- intersect(L4AlphaEGO@result$Description, L4LambdaEGO@result$Description )                # identify common genes
L4_alpha <- L4AlphaEGO@result$Description[!(L4AlphaEGO@result$Description %in% L4LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L4_lambda <- L4LambdaEGO@result$Description[!(L4LambdaEGO@result$Description %in% L4AlphaEGO@result$Description)]
L4_lambda <- L4LambdaEGO[L4LambdaEGO@result$Description %in% L4_lambda, ] 
L4_lambda_list <- L4_lambda[, c("ID", "p.adjust")] 
write_tsv(L4_lambda_list, "honours/results/DEAnalysis/FortopGO/L4_lambda_list.tsv")

L3_common <- intersect(L3AlphaEGO@result$Description, L3LambdaEGO@result$Description )                # identify common genes
L3_alpha <- L3AlphaEGO@result$Description[!(L3AlphaEGO@result$Description %in% L3LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L3_lambda <- L3LambdaEGO@result$Description[!(L3LambdaEGO@result$Description %in% L3AlphaEGO@result$Description)]
L3_lambda <- L3LambdaEGO[L3LambdaEGO@result$Description %in% L3_lambda, ] 
L3_lambda_list <- L3_lambda[, c("ID", "p.adjust")] 
write_tsv(L3_lambda_list, "honours/results/DEAnalysis/FortopGO/L3_lambda_list.tsv")

L2_common <- intersect(L2AlphaEGO@result$Description, L2LambdaEGO@result$Description )                # identify common genes
L2_alpha <- L2AlphaEGO@result$Description[!(L2AlphaEGO@result$Description %in% L2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L2_lambda <- L2LambdaEGO@result$Description[!(L2LambdaEGO@result$Description %in% L2AlphaEGO@result$Description)]
L2_lambda <- L2LambdaEGO[L2LambdaEGO@result$Description %in% L2_lambda, ] 
L2_lambda_list <- L2_lambda[, c("ID", "p.adjust")] 
write_tsv(L2_lambda_list, "honours/results/DEAnalysis/FortopGO/L2_lambda_list.tsv")

L1_common <- intersect(L1AlphaEGO@result$Description, L1LambdaEGO@result$Description )                # identify common genes
L1_alpha <- L1AlphaEGO@result$Description[!(L1AlphaEGO@result$Description %in% L1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
L1_lambda <- L1LambdaEGO@result$Description[!(L1LambdaEGO@result$Description %in% L1AlphaEGO@result$Description)]
L1_lambda <- L1LambdaEGO[L1LambdaEGO@result$Description %in% L1_lambda, ] 
L1_lambda_list <- L1_lambda[, c("ID", "p.adjust")] 
write_tsv(L1_lambda_list, "honours/results/DEAnalysis/FortopGO/L1_lambda_list.tsv")

M1_common <- intersect(M1AlphaEGO@result$Description, M1LambdaEGO@result$Description )                # identify common genes
M1_alpha <- M1AlphaEGO@result$Description[!(M1AlphaEGO@result$Description %in% M1LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
M1_lambda <- M1LambdaEGO@result$Description[!(M1LambdaEGO@result$Description %in% M1AlphaEGO@result$Description)]
M1_lambda <- M1LambdaEGO[M1LambdaEGO@result$Description %in% M1_lambda, ] 
M1_lambda_list <- M1_lambda[, c("ID", "p.adjust")] 
write_tsv(M1_lambda_list, "honours/results/DEAnalysis/FortopGO/M1_lambda_list.tsv")

M2_common <- intersect(M2AlphaEGO@result$Description, M2LambdaEGO@result$Description )                # identify common genes
M2_alpha <- M2AlphaEGO@result$Description[!(M2AlphaEGO@result$Description %in% M2LambdaEGO@result$Description)] # - (all alpha in lambda) = all alpha not in lambda
M2_lambda <- M2LambdaEGO@result$Description[!(M2LambdaEGO@result$Description %in% M2AlphaEGO@result$Description)]
M2_lambda <- M2LambdaEGO[M2LambdaEGO@result$Description %in% M2_lambda, ] 
M2_lambda_list <- M2_lambda[, c("ID", "p.adjust")] 
write_tsv(M2_lambda_list, "honours/results/DEAnalysis/FortopGO/M2_lambda_list.tsv")

length(M1_alpha)
dim(M1_lambda)

# KEGG pathway analysis 





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

