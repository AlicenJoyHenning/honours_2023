# Differential Expression analysis in R (seurat

###### dependencies & load object #####
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

treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")

TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'Mono',
                                   '1' = 'Mono',
                                   '2' = 'CD4_helper',
                                   '3' = 'CD4_naive',
                                   '4' = 'Neutro',
                                   '5' = 'CD8',
                                   '6' = 'B',
                                   '7' = 'CD8',
                                   '8' = 'CD8',
                                   '9' = 'NK',
                                   '10' = '10',
                                   '11' = 'Tregs',
                                   '12' = '12',
                                   '13' = '13',
                                   '14' = '14')



##### View the annotated clusters #####
###
Idents(TreatmentAnnotated)
palette.b <- c(
               "#d72554", #0 and 1
               "#6ab5ba", #2
               "#2e8f95", #3
               "#900c3e", #4
               "#8caf2e", #5
               "#c35cad", #6
               "#3b7749", #9
               "#FFCB3E", #10
               "#00a68e", #11
               "#5c040c", #12
               "#f7bc6e", #13
               "#6ab5ba" #13 
               )
all <- DimPlot(TreatmentAnnotated, reduction = "umap", pt.size = 1.5, label = TRUE, label.color = "white", label.size = 6, label.box = TRUE, repel = TRUE, cols = palette.b)

features <- c("CXCR1", "CXCR2", "CD4", "KLF2", "IL7R", "CD8B", "CD8A", "BCL11A", "CD40", "NKG7", "CTSW", "FOXP3", "RTKN2")
f <- DotPlot(object = treatment, 
        features = features,
        cols = c("grey", "#6ab5ba")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

all | f

##### DE ######

# {1} : Creating metadata column to hold cell type AND stimulation information (as opposed to just cell type) : 

TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$stim, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the idents to that column 

# Now we can easily refer to a cell & treatment type, to see the options : 

levels(TreatmentAnnotated) # list the options to perform DE on :

# MYELOID CELL TYPES : 
# 1. "Mono_alpha", "Mono_lambda"           VS    Mono_untreated
# 2. "Neutro_alpha", Neutro_lambda"        VS    Neutro_untreated 

# LYMPHOID CELL TYPES: 
# T CELLS : 
# 1. "CD4_helper_alpha", "CD4_helper_lambda" VS    "CD4_helper_untreated" 
# 2. "CD4_naive_alpha", "CD4_naive_lambda"   VS    "CD4_naive_untreated" 
# 3. "Tregs_alpha", "Tregs_lambda"           VS    "Tregs_untreated"
# 4. "CD8_alpha", "CD8_lambda"               VS    "CD8_untreated"
# 5. "NK_alpha", "NK_lambda"                 VS    "NK_untreated" 

# B CELLS : 
# 6. "B_alpha", B_lambda"                    VS    "B_untreated"      
    

# {2} Use FindMarkers() to find the genes that are different between stimulated and untreated cell types
?FindMarkers()


# MYELOID CELL TYPES : 

# 1. Monocytes 
M1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Mono_alpha", ident.2 = "Mono_untreated", sep = "_") 
M1AlphaResponse$gene <- rownames(M1AlphaResponse) # puts gene names into a column
M1AlphaResponse <- filter(M1AlphaResponse, p_val_adj < 0.05)
down <- filter(M1AlphaResponse, avg_log2FC < 0)
up <- filter(M1AlphaResponse, avg_log2FC >= 0)
MonoAlphaResponse <- data.frame(Gene = M1AlphaResponse$gene, Log2FoldChange = M1AlphaResponse$avg_log2FC)  # for cluster profiler 
write.csv(MonoAlphaResponse, "honours/results/DEAnalysis/after_filtering/MonoAlphaResponse.csv", row.names = FALSE)


M1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Mono_lambda", ident.2 = "Mono_untreated", sep = "_") 
M1LambdaResponse$gene <- rownames(M1LambdaResponse) 
M1LambdaResponse <- filter(M1LambdaResponse, p_val_adj < 0.05)
down <- filter(M1LambdaResponse, avg_log2FC < 0)
up <- filter(M1LambdaResponse, avg_log2FC >= 0)
MonoLambdaResponse <- data.frame(Gene = M1LambdaResponse$gene, Log2FoldChange = M1LambdaResponse$avg_log2FC)  
write.csv(MonoLambdaResponse, "honours/results/DEAnalysis/after_filtering/MonoLambdaResponse.csv", row.names = FALSE)


# 2. Neutrophils 

M2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Neutro_alpha", ident.2 = "Neutro_untreated", sep = "_") 
M2AlphaResponse$gene <- rownames(M2AlphaResponse)
M2AlphaResponse <- filter(M2AlphaResponse, p_val_adj < 0.05)
down <- filter(M2AlphaResponse, avg_log2FC < 0)
up <- filter(M2AlphaResponse, avg_log2FC >= 0)
neuAlphaResponse <- data.frame(Gene = M2AlphaResponse$gene, Log2FoldChange = M2AlphaResponse$avg_log2FC)  
write.csv(neuAlphaResponse, "honours/results/DEAnalysis/after_filtering/neuAlphaResponse.csv", row.names = FALSE)
# size = 1507 DEGs 

M2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Neutro_lambda", ident.2 = "Neutro_untreated", sep = "_") 
M2LambdaResponse$gene <- rownames(M2LambdaResponse) 
M2LambdaResponse <- filter(M2LambdaResponse, p_val_adj < 0.05)
down <- filter(M2LambdaResponse, avg_log2FC < 0)
up <- filter(M2LambdaResponse, avg_log2FC >= 0)
neuLambdaResponse <- data.frame(Gene = M2LambdaResponse$gene, Log2FoldChange = M2LambdaResponse$avg_log2FC)  
write.csv(neuLambdaResponse, "honours/results/DEAnalysis/after_filtering/neuLambdaResponse.csv", row.names = FALSE)


# LYMPHOID CELL TYPES : 
# 1. CD4_helper T cells 
L1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_helper_alpha", ident.2 = "CD4_helper_untreated", sep = "_") 
L1AlphaResponse$gene <- rownames(L1AlphaResponse) 
L1AlphaResponse  <- filter(L1AlphaResponse , p_val_adj < 0.05)
down <- filter(L1AlphaResponse, avg_log2FC < 0)
up <- filter(L1AlphaResponse, avg_log2FC >= 0)
CD4hAlphaResponse <- data.frame(Gene = L1AlphaResponse$gene, Log2FoldChange = L1AlphaResponse$avg_log2FC)  
write.csv(CD4hAlphaResponse, "honours/results/DEAnalysis/after_filtering/CD4hAlphaResponse.csv", row.names = FALSE)


L1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_helper_lambda", ident.2 = "CD4_helper_untreated", sep = "_") 
L1LambdaResponse$gene <- rownames(L1LambdaResponse)
L1LambdaResponse <- filter(L1LambdaResponse, p_val_adj < 0.05)
down <- filter(L1LambdaResponse, avg_log2FC < 0)
up <- filter(L1LambdaResponse, avg_log2FC >= 0)
CD4hLambdaResponse <- data.frame(Gene = L1LambdaResponse$gene, Log2FoldChange = L1LambdaResponse$avg_log2FC)  
write.csv(CD4hLambdaResponse, "honours/results/DEAnalysis/after_filtering/CD4hLambdaResponse.csv", row.names = FALSE)


# 2. CD4_naive T cells 

L2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_naive_alpha", ident.2 = "CD4_naive_untreated", sep = "_") 
L2AlphaResponse$gene <- rownames(L2AlphaResponse) 
L2AlphaResponse <- filter(L2AlphaResponse, p_val_adj < 0.05)
down <- filter(L2AlphaResponse, avg_log2FC < 0)
up <- filter(L2AlphaResponse, avg_log2FC >= 0)
CD4nAlphaResponse <- data.frame(Gene = L2AlphaResponse$gene, Log2FoldChange = L2AlphaResponse$avg_log2FC)  
write.csv(CD4nAlphaResponse, "honours/results/DEAnalysis/after_filtering/CD4nAlphaResponse.csv", row.names = FALSE)


L2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_naive_lambda", ident.2 = "CD4_naive_untreated", sep = "_") 
L2LambdaResponse$gene <- rownames(L2LambdaResponse)
L2LambdaResponse <- filter(L2LambdaResponse, p_val_adj < 0.05)
down <- filter(L2LambdaResponse, avg_log2FC < 0)
up <- filter(L2LambdaResponse, avg_log2FC >= 0)
CD4nLambdaResponse <- data.frame(Gene = L2LambdaResponse$gene, Log2FoldChange = L2LambdaResponse$avg_log2FC)  
write.csv(CD4nLambdaResponse, "honours/results/DEAnalysis/after_filtering/CD4nLambdaResponse.csv", row.names = FALSE)
 

# 3. Tregs_ 
L3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_alpha", ident.2 = "Tregs_untreated", sep = "_") 
L3AlphaResponse$gene <- rownames(L3AlphaResponse)
L3AlphaResponse <- filter(L3AlphaResponse, p_val_adj < 0.05)
down <- filter(L3AlphaResponse, avg_log2FC < 0)
up <- filter(L3AlphaResponse, avg_log2FC >= 0)
TregsAlphaResponse <- data.frame(Gene = L3AlphaResponse$gene, Log2FoldChange = L3AlphaResponse$avg_log2FC)  
write.csv(TregsAlphaResponse, "honours/results/DEAnalysis/after_filtering/TregsAlphaResponse.csv", row.names = FALSE)
 

L3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_lambda", ident.2 = "Tregs_untreated", sep = "_") 
L3LambdaResponse$gene <- rownames(L3LambdaResponse)
L3LambdaResponse <- filter(L3LambdaResponse, p_val_adj < 0.05)
down <- filter(L3LambdaResponse, avg_log2FC < 0)
up <- filter(L3LambdaResponse, avg_log2FC >= 0)
TregsLambdaResponse <- data.frame(Gene = L3LambdaResponse$gene, Log2FoldChange = L3LambdaResponse$avg_log2FC)  
write.csv(TregsLambdaResponse, "honours/results/DEAnalysis/after_filtering/TregsLambdaResponse.csv", row.names = FALSE)


# 4. CD8_alpha 
L4AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD8_alpha", ident.2 = "CD8_untreated", sep = "_") 
L4AlphaResponse$gene <- rownames(L4AlphaResponse)
L4AlphaResponse <- filter(L4AlphaResponse, p_val_adj < 0.05)
down <- filter(L4AlphaResponse, avg_log2FC < 0)
up <- filter(L4AlphaResponse, avg_log2FC >= 0)
CD8AlphaResponse <- data.frame(Gene = L4AlphaResponse$gene, Log2FoldChange = L4AlphaResponse$avg_log2FC)  
write.csv(CD8AlphaResponse, "honours/results/DEAnalysis/after_filtering/CD8AlphaResponse.csv", row.names = FALSE)
 

L4LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD8_lambda", ident.2 = "CD8_untreated", sep = "_") 
L4LambdaResponse$gene <- rownames(L4LambdaResponse) 
L4LambdaResponse <- filter(L4LambdaResponse, p_val_adj < 0.05)
down <- filter(L4LambdaResponse, avg_log2FC < 0)
up <- filter(L4LambdaResponse, avg_log2FC >= 0)
CD8LambdaResponse <- data.frame(Gene = L4LambdaResponse$gene, Log2FoldChange = L4LambdaResponse$avg_log2FC)  
write.csv(CD8LambdaResponse, "honours/results/DEAnalysis/after_filtering/CD8LambdaResponse.csv", row.names = FALSE)
 

# 5. NK cells 
L5AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_alpha", ident.2 = "NK_untreated", sep = "_") 
L5AlphaResponse$gene <- rownames(L5AlphaResponse)
L5AlphaResponse <- filter(L5AlphaResponse, p_val_adj < 0.05)
down <- filter(L5AlphaResponse, avg_log2FC < 0)
up <- filter(L5AlphaResponse, avg_log2FC >= 0)
NKAlphaResponse <- data.frame(Gene = L5AlphaResponse$gene, Log2FoldChange = L5AlphaResponse$avg_log2FC)  
write.csv(NKAlphaResponse, "honours/results/DEAnalysis/after_filtering/NKAlphaResponse.csv", row.names = FALSE)


L5LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_lambda", ident.2 = "NK_untreated", sep = "_") 
L5LambdaResponse$gene <- rownames(L5LambdaResponse)
L5LambdaResponse <- filter(L5LambdaResponse, p_val_adj < 0.05)
down <- filter(L5LambdaResponse, avg_log2FC < 0)
up <- filter(L5LambdaResponse, avg_log2FC >= 0)
NKLambdaResponse <- data.frame(Gene = L5LambdaResponse$gene, Log2FoldChange = L5LambdaResponse$avg_log2FC)  
write.csv(CD8LambdaResponse, "honours/results/DEAnalysis/after_filtering/NKLambdaResponse.csv", row.names = FALSE)



# 6. B CELLS 
L6AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_alpha", ident.2 = "B_untreated", sep = "_") # size = 1287
L6AlphaResponse$gene <- rownames(L6AlphaResponse) # puts gene names into a column
L6AlphaResponse <- filter(L6AlphaResponse, p_val_adj < 0.05)
down <- filter(L6AlphaResponse, avg_log2FC < 0)
up <- filter(L6AlphaResponse, avg_log2FC >= 0)
BAlphaResponse <- data.frame(Gene = L6AlphaResponse$gene, Log2FoldChange = L6AlphaResponse$avg_log2FC)  
write.csv(BAlphaResponse, "honours/results/DEAnalysis/after_filtering/BAlphaResponse.csv", row.names = FALSE) # next step requires entry as a csv hence this step 


L6LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_lambda", ident.2 = "B_untreated", sep = "_") # size = 310 
L6LambdaResponse$gene <- rownames(L6LambdaResponse)
L6LambdaResponse <- filter(L6LambdaResponse, p_val_adj < 0.05)
down <- filter(L6LambdaResponse, avg_log2FC < 0)
up <- filter(L6LambdaResponse, avg_log2FC >= 0)
BLambdaResponse <- data.frame(Gene = L6LambdaResponse$gene,Log2FoldChange = L6LambdaResponse$avg_log2FC)  
write.csv(BLambdaResponse, "honours/results/DEAnalysis/after_filtering/BLambdaResponse.csv", row.names = FALSE)


##### 
##### Extra use of FindMarkers() ? #####

# [a] To find common DEG between 2 treatments compared to the control (same pathways?) : 
B_cells_common_genes <- intersect(AlphaResponse$gene, LambdaResponse$gene)                # identify common genes
B_cells_common_alpha <- AlphaResponse[AlphaResponse$gene %in% B_cells_common_genes, ]     # Extract common genes from each dataframe 
B_cells_common_lambda <- LambdaResponse[LambdaResponse$gene %in% B_cells_common_genes, ]
B_cells_common_dataframe <- merge(B_cells_common_alpha, B_cells_common_lambda, by = 'gene', all = TRUE) # # merge together the 2 subsetted dataframes to create 1 containing only common genes :
# size = 223
# [b] To find treatment specific DEGs : 
B_cells_alpha_genes <- AlphaResponse$gene[!(AlphaResponse$gene %in% LambdaResponse$gene)] # - (all alpha in lambda) = all alpha not in lambda
B_cells_alpha <- AlphaResponse[AlphaResponse$gene %in% B_cells_alpha_genes, ]
# size = 1064 genes unique to alpha 
B_cells_lambda_genes <- LambdaResponse$gene[!(LambdaResponse$gene %in% AlphaResponse$gene)]
B_cells_lambda <- LambdaResponse[LambdaResponse$gene %in% B_cells_lambda_genes, ] 
# size = 87 genes unique to lambda
# [c] Store information in excel file 
B_cells <- list('B_cells_Alpha_Response' = AlphaResponse, 'B_cells_Lambda_Response' = LambdaResponse,'B_cells_common_response' = B_cells_common_dataframe, 'B_cells_alpha_specific' = B_cells_alpha,'B_cells_lambda_specific' = B_cells_lambda)
openxlsx::write.xlsx(B_cells, file = "honours/results/DEAnalysis/B_cells.xlsx")

##### Visualisations ##### 
# TIER ONE : general trends (no stats)

# [a] Histogram showing DEGs across cell types and treatments (9 sections)

##
B_FP <- FeaturePlot(TreatmentAnnotated, 
features = c("ISG15"),# "IFIT3",# "IFI6",# "IFIT2",# "MX1",# "TNFSF10",# "SAMD9L",# "HERC5",# "IFIT1",# "IFI44L"), ),
split.by = "stim",            
cols = c("grey", "red")) + geom_hline(yintercept = c(-14, -11), linetype = "dotted", color = "black") 


theme_set(theme_cowplot())
t.cells <- subset(TreatmentAnnotated, idents = "Naive CD4+ T cells")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)
p1 <- ggplot(avg.t.cells, 
             aes(untreated, alpha)) + 
  geom_point() +
  xlab("UNTREATED") +
  ylab("ALPHA")

p2 <- ggplot(avg.t.cells, 
             aes(untreated, lambda)) +
  geom_point() +
  xlab("UNTREATED") +
  ylab("LAMBDA")
p3 <- p1 + p2 
p4 <- p3 +  ggtitle("Differentially expressed genes in Naive T cells") 
p4

##### ClusterProfiler #####

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





##### KEGG analysis #####

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




##### DGEA Plots #####
library(enrichplot)

# Plot the enrichment results
barplot(L6_lambda$qvalue, showCategory=20)



###### pplots not yet uesd #####

egox <- setReadable(topego, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(egox, FoldChange = de)
p3 <- cnetplot(egox, FoldChange= de, circular = TRUE, colorEdge = TRUE) 

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox)
p2 <- cnetplot(ego, foldChange=gene_list)
p3 <- cnetplot(edox, foldChange=gene_list, circular = TRUE, colorEdge = TRUE) 
p4 <- heatplot(edox, foldChange=gene_list, showCategory=5)

edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# 2 : Prepare Input : read in DE csv file 
BCellAlphadf = read.csv("honours/results/DEAnalysis/BCellAlphaResponse.csv", header = TRUE)

gene_list <- BCellAlphadf$Log2FoldChange
names(gene_list) <- BCellAlphadf$Gene
gene_list = sort(gene_list, decreasing = TRUE)

# 3 : Gene Set ENrichment 

keytypes(org.Hs.eg.db) # see what options are available 

gse <- gseGO( geneList = gene_list, 
              ont = "ALL", # enrichment on all 3 GO branches (biological processes, molecular function, cellular component)
              keyType = "SYMBOL", 
              nPerm = 10000, 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = organism, 
              pAdjustMethod = "none")
goplot(gse)


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

##### EnrichR ######
BiocManager::install("enrichR")
library("enrichR")
library(enrichplot)
?DEenrichRPlot()

DEenrichRPlot(object = TreatmentAnnotated,
  ident.1 = "B_cells_alpha",
  ident.2 = "B_cells_untreated",
  balanced = TRUE,
  logfc.threshold = 0.25,
  assay = "RNA",
  max.genes = 10000,
  test.use = "wilcox",
  p.val.cutoff = 0.05,
  cols = NULL,
  enrich.database = "GO_All",
  num.pathway = 20,
  return.gene.list = FALSE,
)

"#6ab5ba""#6ab5ba"
##### TopGO #####
BiocManager::install("topGO")
library(topGO)



