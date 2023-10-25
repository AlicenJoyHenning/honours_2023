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
library(dplyr) # data frame work 
library(patchwork)
library(openxlsx)
install.packages("dbplyr")
library(dbplyr)
library(readr)

library(ggplot2)
library(reshape2)

# [1] Alter Seurat object ####
TreatmentAnnotated <- read_rds("honours/results/FinalIndex/TreatmentAnnotated.rds")
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
getwd()

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

  # Create R object suitable for analysis : (for merging)
  EGOdf <- data.frame(
    Type = rep(paste0(CellType), nrow(EGO)),
    Treatment = rep(paste0(Treatment), nrow(EGO)),
    Regulation = rep(paste0(Regulation), nrow(EGO)),
    GODescription = EGO@result$Description,
    GeneRatio = EGO@result$GeneRatio
  ) 
  saveRDS(EGOdf, file = paste0("GO_", CellType, "_", Treatment, "_", Regulation, ".rds")) 
  cat("Downstream Dataframe created", "\n")
  
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

# [4] Identify and quantify GO terms ####

n.a.u <- readRDS("GO_neutrophils_alpha_upregulated.rds") 
n.a.d <- readRDS("GO_neutrophils_alpha_downregulated.rds")
n.l.u <- readRDS("GO_neutrophils_lambda_upregulated.rds")
n.l.d <- readRDS("GO_neutrophils_lambda_downregulated.rds")

m.a.u <- readRDS("GO_myeloid_alpha_upregulated.rds") 
m.a.d <- readRDS("GO_myeloid_alpha_downregulated.rds")
m.l.u <- readRDS("GO_myeloid_lambda_upregulated.rds")
m.l.d <- readRDS("GO_myeloid_lambda_downregulated.rds")

b.a.u <- readRDS("GO_B_alpha_upregulated.rds") 
b.a.d <- readRDS("GO_B_alpha_downregulated.rds")
b.l.u <- readRDS("GO_B_lambda_upregulated.rds")
b.l.d <- readRDS("GO_B_lambda_downregulated.rds")

l.a.u <- readRDS("GO_lymphoid_alpha_upregulated.rds") 
l.a.d <- readRDS("GO_lymphoid_alpha_downregulated.rds")
l.l.u <- readRDS("GO_lymphoid_lambda_upregulated.rds")
l.l.d <- readRDS("GO_lymphoid_lambda_downregulated.rds")

GOanalysis <-  bind_rows(n.a.u, n.a.d, n.l.u, n.l.d, 
                         m.a.u, m.a.d, m.l.u, m.l.d, 
                         b.a.u, b.a.d, b.l.u, b.l.d, 
                         l.a.u, l.a.d, l.l.u, l.l.d)   

GODownAnalysis <- bind_rows(n.a.d, n.l.d, 
                            m.a.d, m.l.d, 
                            b.a.d, b.l.d, 
                            l.a.d, l.l.d)


# [4] Heat Map itself ####

# UPREGULATED

# Transpose the Data Frame
heatT <- data.frame(
  GO = c("GO1", "GO2", "GO3", "GO4", "GO5"),
  neutrophils_α = c(40, 36, 44, 19, 30),
  neutrophils_λ = c(0,  0,  0,  0,  0),
  myeloid_α =     c(46, 35, 48, 47, 38),
  myeloid_λ =     c(0,  0,  0,  0,  0),
  Bcells_α =      c(38, 37, 40, 13, 28),
  Bcells_λ =      c(20, 18, 14,  0, 12),
  lymphoid_α =    c(42, 34, 46, 26, 35),
  lymphoid_λ =    c(5,  3,  0,  0,  0)
)

heatTalt <- data.frame(
  GO = c("GO1", "GO2", "GO3", "GO4", "GO5"),
  neutrophils_α = c(40, 36, 44, 19, 30),
  myeloid_α =     c(46, 35, 48, 47, 38),
  Bcells_α =      c(38, 37, 40, 13, 28),
  lymphoid_α =    c(42, 34, 46, 26, 35),
  neutrophils_λ = c(0,  0,  0,  0,  0),
  myeloid_λ =     c(0,  0,  0,  0,  0),
  Bcells_λ =      c(20, 18, 14,  0, 12),
  lymphoid_λ =    c(5,  3,  0,  0,  0)
)


# Reorder the GO column 
heatTalt$GO <- factor(heatTalt$GO, levels = c("GO6", "GO5", "GO4", "GO3", "GO2", "GO1"))


heatTalt_long <- melt(heatTalt, id.vars = "GO", variable.name = "CellType", value.name = "Count")


up <- ggplot(heatTalt_long, aes(x = CellType, y = GO, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgrey", high = "#6ab5ba") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold") 
        ) +  # Rotate x-axis labels and set font size and style
  labs(x = "", y = "Top enriched GO Terms", fill = "Count") +
  theme(legend.position = "left") +  
  # scale_x_discrete(labels = c("Bcells_α" = "B Cells α", 
  #                             "Bcells_λ" = "B Cells λ", 
  #                             "lymphoid_α" = "Lymphoid α", 
  #                             "lymphoid_λ" = "Lymphoid λ",
  #                             "neutrophils_α" = "Neutrophils α",
  #                             "neutrophils_λ" = "Neutrophils λ",
  #                             "myeloid_α" = "Myeloid α",
  #                             "myeloid_λ" = "Myeloid λ")) +  
  scale_y_discrete(labels = c("GO1" = "1", 
                              "GO2" = "2", 
                              "GO3" = "3", 
                              "GO4" = "4",
                              "GO5" = "5"))

# DOWN REGULATED
# Transpose the Data Frame
heatT <- data.frame(
  GO = c("1", "2", "3", "4", "5"),
  neutrophils_α = c(0, 25, 20, 12, 14),
  neutrophils_λ = c(1,  0,  0,  0,  0),
  myeloid_α =     c(0, 31, 18,16, 10),
  myeloid_λ =     c(1, 0,  0,  0,  0),
  Bcells_α =      c(18, 10, 0, 0, 0),
  Bcells_λ =      c(1, 0, 0, 0,  0),
  lymphoid_α =    c(43, 12, 0, 0, 0),
  lymphoid_λ =    c(1, 0,  0,  0,  0)
)

heatTaltt <- data.frame(
  GO = c("1", "2", "3", "4", "5"),
  neutrophils_α = c(0, 25, 20, 12, 14),
  myeloid_α =     c(0, 31, 18,16, 10),
  Bcells_α =      c(18, 10, 0, 0, 0),
  lymphoid_α =    c(43, 12, 0, 0, 0),
  
  
  neutrophils_λ = c(1,  0,  0,  0,  0),
  myeloid_λ =     c(1, 0,  0,  0,  0),
  Bcells_λ =      c(1, 0, 0, 0,  0),
  lymphoid_λ =    c(1, 0,  0,  0,  0)
)

# Reorder the GO column 
heatTaltt$GO <- factor(heatTaltt$GO, levels = c("5", "4", "3", "2", "1"))


heatTaltt_long <- melt(heatTaltt, id.vars = "GO", variable.name = "CellType", value.name = "Count")


down <- ggplot(heatTaltt_long, aes(x = CellType, y = GO, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgrey", high = "#868686") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)  #, face = "bold") 
  ) +  # Rotate x-axis labels and set font size and style
  labs(x = "", y = "Top enriched GO Terms", fill = "Counts") +
  theme(legend.position = "left") +  
  scale_x_discrete(labels = c("Bcells_α" = "B Cells α", 
                              "Bcells_λ" = "B Cells λ", 
                              "lymphoid_α" = "Lymphoid α", 
                              "lymphoid_λ" = "Lymphoid λ",
                              "neutrophils_α" = "Neutrophils α",
                              "neutrophils_λ" = "Neutrophils λ",
                              "myeloid_α" = "Myeloid α",
                              "myeloid_λ" = "Myeloid λ")) +  
  scale_y_discrete(labels = c("GO1" = "1", 
                              "GO2" = "2", 
                              "GO3" = "3", 
                              "GO4" = "4",
                              "GO5" = "5")) 







# down with altered positioning (aesthetic)
down <- ggplot(heatTaltt_long, aes(x = CellType, y = GO, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightgrey", high = "#868686") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)  #, face = "bold") 
  ) +  # Rotate x-axis labels and set font size and style
  labs(x = "", y = "Top enriched GO Terms", fill = "Counts") +
  theme(legend.position = "left") +  
  scale_x_discrete(labels = c("Bcells_α" = "B Cells", 
                              "Bcells_λ" = "B Cells", 
                              "lymphoid_α" = "Lymphoid", 
                              "lymphoid_λ" = "Lymphoid",
                              "neutrophils_α" = "Neutrophils",
                              "neutrophils_λ" = "Neutrophils",
                              "myeloid_α" = "Myeloid",
                              "myeloid_λ" = "Myeloid")) +  
  scale_y_discrete(labels = c("GO1" = "1", 
                              "GO2" = "2", 
                              "GO3" = "3", 
                              "GO4" = "4",
                              "GO5" = "5"))
down
up / down
