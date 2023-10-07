# create heat map showing the top up and down regulated GO terms for broad cell populations 
# 1) neutrophils, 2) myeloid, 3) B cells, 4) lymphoid

# [1] Alter Seurat object ####
SubTreatment <- subset(TreatmentAnnotated, seurat_clusters != 13 & seurat_clusters != 17 & seurat_clusters != 12)

SubTreatment <- RenameIdents(SubTreatment, 
                                   'monocytes' = 'myeloid',
                                   'naive CD4 T' = 'lymphoid',
                                   'T helper' = 'lymphoid',
                                   'naive CD8 T'= 'lymphoid',
                                   'cytotoxic T'= 'lymphoid',
                                   'mDCs' = 'myeloid',
                                   'NKT'= 'lymphoid',
                                   'Tregs'= 'lymphoid',
                                   'NK' = 'lymphoid',
                                   'DCs' = 'myeloid', 
                                   'Tcm' = 'lymphoid')

SubTreatment$celltype.stim <- paste(Idents(SubTreatment), SubTreatment$sample, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
SubTreatment$celltype <- Idents(SubTreatment) # restores the cell type column 
Idents(SubTreatment) <- "celltype.stim"

setwd("honours/mydreamheatmap/")

# [2] DEG on categories above for alpha and lambda treatment  ####


# FUNCTION

DGEA <- function(CellType) {
  
  # ALPHA
  # find differentially expressed genes :
  ACellType <- FindMarkers(SubTreatment, ident.1 = paste0(CellType, "_alpha"), ident.2 = paste0(CellType, "_untreated"), sep = "_")
  ACellType$gene <- rownames(ACellType) # puts gene names into a column
  
  # find significant genes :
  ACellType <- filter(ACellType, p_val_adj < 0.05)
  Atotal <- nrow(ACellType)
  
  # identify which of these are up and down regulated :
  CellTypeAdown <- filter(ACellType, avg_log2FC < 0)
  Adown <- nrow(CellTypeAdown)
  CellTypeAup <- filter(ACellType, avg_log2FC >= 0)
  Aup <- nrow(CellTypeAup)
  
  # create text file outputs with this info :
  write.table(CellTypeAdown$gene, file = paste0(CellType, "_alpha_downregulated.txt"), row.names = FALSE, col.names = FALSE)
  write.table(CellTypeAup$gene, file = paste0(CellType, "_alpha_upregulated.txt"), row.names = FALSE, col.names = FALSE)
  
  
  
  # LAMBDA
  # find differentially expressed genes :
  LCellType <- FindMarkers(SubTreatment, ident.1 = paste0(CellType, "_lambda"), ident.2 = paste0(CellType, "_untreated"), sep = "_")
  LCellType$gene <- rownames(LCellType) # puts gene names into a column
  
  # find significant genes :
  LCellType <- filter(LCellType, p_val_adj < 0.05)
  Ltotal <- nrow(LCellType)
  
  # identify which of these are up and down regulated :
  CellTypeLdown <- filter(LCellType, avg_log2FC < 0)
  Ldown <- nrow(CellTypeLdown)
  CellTypeLup <- filter(LCellType, avg_log2FC >= 0)
  Lup <- nrow(CellTypeLup)
  
  # create text file outputs with this info :
  write.table(CellTypeLdown$gene, file = paste0( CellType,"_lambda_downregulated.txt"), row.names = FALSE, col.names = FALSE)
  write.table(CellTypeLup$gene, file = paste0(CellType, "_lambda_upregulated.txt"), row.names = FALSE, col.names = FALSE)
  
  # OUTPUT
  cat("Total DEGs in", CellType, "in response to alpha:", Atotal, "\n")
  cat("UP:", Aup, "\n")
  cat("DOWN:", Adown, "\n")
  
  cat("Total DEGs in", CellType, "in response to lambda:", Ltotal, "\n")
  cat("UP:", Lup, "\n")
  cat("DOWN:", Ldown, "\n")
}

DGEA("lymphoid")


# [3] GO analysis on generated DEG lists (per category and up and down) ####

# [4] heatmap itself ####