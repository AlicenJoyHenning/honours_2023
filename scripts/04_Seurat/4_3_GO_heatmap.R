# create heat map showing the top up and down regulated GO terms for broad cell populations 
# 1) neutrophils, 2) myeloid, 3) B cells, 4) lymphoid

# [0] Dependencies ####

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
install.packages("dbplyr")
library(dbplyr)
library(readr)

library(ggplot2)
library(reshape2)

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
  saveRDS(CellTypeAdown, paste0(CellType,"_alpha_downregulated.rds"))
  saveRDS(CellTypeAup, paste0(CellType,"_alpha_upregulated.rds"))
  
  
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
  saveRDS(CellTypeLdown, paste0(CellType,"_lambda_downregulated.rds"))
  saveRDS(CellTypeLup, paste0(CellType,"_lambda_upregulated.rds"))
  
  # OUTPUT
  cat("Total DEGs in", CellType, "in response to alpha:", Atotal, "\n")
  cat("UP:", Aup, "\n")
  cat("DOWN:", Adown, "\n")
  
  cat("Total DEGs in", CellType, "in response to lambda:", Ltotal, "\n")
  cat("UP:", Lup, "\n")
  cat("DOWN:", Ldown, "\n")
}

DGEA("neutrophils")
DGEA("myeloid")
DGEA("B")
DGEA("lymphoid")

# [3] GO analysis on generated DEG lists (per category and up and down) ####

# 3.1 : install and load annotation for desired organism (human) : org.Hs.eg.db
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# 3.2 : have to convert DE gene list to contain entrezgene_ID BIOMART
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Create a Mart object for the human genome
# listDatasets(mart) # List available datasets for the Mart (214!)

# 3.3 :  Retrieve gene information for human genes :
gene_ID <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
                 mart = mart) 



# 3.4 Function 
EvenMoreEnrichedGO <- function(CellType, Treatment, Regulation) {
  # read in DEGs R object : 
  DEGs <- readRDS(paste0("DEGs/", CellType, "_", Treatment, "_", Regulation, ".rds"))
  cat("DEGs read in", "\n")
  # merge with mart object to get entrez gene names : 
  DEGs <- unique(merge(DEGs, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name"))
  cat("DEGs merged", "\n")
  # create object to perform GO analysis on, just the entrezgene  id column :
  DEGs <- na.omit(DEGs$entrezgene_id)
  
  # perform GO analysis :
  EGO <- enrichGO(gene = DEGs, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",ont = "BP",pAdjustMethod = "BH", pvalueCutoff = 0.05)
  EGO <- filter(EGO, EGO@result$p.adjust < 0.05)
  cat("enrich GO performed", "\n")
  # save as R object and create text file with GO term list : 
  output <- paste(EGO$Description, EGO$GeneRatio, sep = "\t") # concatenate values from 2 columns 
  write.table(output, file = paste0("GO_", CellType, "_", Treatment, "_", Regulation, ".txt"), row.names = FALSE, col.names = FALSE)
  cat("Output saved")
  } 

# neutrophils
EvenMoreEnrichedGO("neutrophils", "alpha", "upregulated")
EvenMoreEnrichedGO("neutrophils", "alpha", "downregulated")
EvenMoreEnrichedGO("neutrophils", "lambda", "upregulated")
EvenMoreEnrichedGO("neutrophils", "lambda", "downregulated")

# myeloid 
EvenMoreEnrichedGO("myeloid", "alpha", "upregulated")
EvenMoreEnrichedGO("myeloid", "alpha", "downregulated")
EvenMoreEnrichedGO("myeloid", "lambda", "upregulated")
EvenMoreEnrichedGO("myeloid", "lambda", "downregulated")

# B cells 
EvenMoreEnrichedGO("B", "alpha", "upregulated")
EvenMoreEnrichedGO("B", "alpha", "downregulated")
EvenMoreEnrichedGO("B", "lambda", "upregulated")
EvenMoreEnrichedGO("B", "lambda", "downregulated")

# lymphoid 
EvenMoreEnrichedGO("lymphoid", "alpha", "upregulated")
EvenMoreEnrichedGO("lymphoid", "alpha", "downregulated")
EvenMoreEnrichedGO("lymphoid", "lambda", "upregulated")
EvenMoreEnrichedGO("lymphoid", "lambda", "downregulated")

# [4] heatmap itself ####

heat <- data.frame(sample = c( 
                         "B cells α", "B cells λ", 
                         "lymphoid α", "lymphoid λ",
                         "neutrophils α", "neutrophils λ",
                         "myeloid α", "myeloid λ"),
                   GO1 = c(30, 20, 17, 4 ,23, 0, 13, 1), 
                   GO2 = c(20, 10, 13, 2, 11, 0, 5, 0), 
                   GO3 = c(16, 10, 12, 0, 15, 0, 10, 0), 
                   GO4 = c(19, 10, 12, 0, 18, 0, 10, 0),
                   GO5 = c(14, 11, 14, 0, 31, 1, 18, 1)) 


heat_long <- melt(heat, id.vars = "sample")

ggplot(heat_long, aes(x = variable, y = sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgrey", high = "#6ab5ba") +
  theme_minimal() +
  labs(x = "GO Terms", y = "Sample", fill = "Count")





# "#8f8e92"