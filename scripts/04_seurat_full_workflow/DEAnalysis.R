# Differential Expression analysis in R (seurat

###### dependencies & load object #####
BiocManager::install("clusterProfiler", version = "3.14")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
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
                                   '0' = 'cMono',
                                   '1' = 'intMono',
                                   '2' = 'CD4_helper',
                                   '3' = 'CD4_naive',
                                   '4' = 'neutrophils',
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
palette.b <- c("#ee5e17", #0
               "#d72554", #1
               "#6ab5ba", #2
               "#2e8f95", #3
               "#900c3e", #4
               "#8caf2e", #5
               "#c35cad", #6
               "#265221", #9
               "#FFCB3E", #10
               "#00a68e", #11
               "#5c040c", #12
               "#f7bc6e", #13
               "#6ab5ba" #14
)
DimPlot(TreatmentAnnotated, reduction = "umap", pt.size = 1.5, label = TRUE, label.color = "white", label.size = 6, label.box = TRUE, repel = TRUE, cols = palette.b)



##### DE ######

# {1} : Creating metadata column to hold cell type AND stimulation information (as opposed to just cell type) : 

TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$stim, sep = "_") # paste puts first entry followed by next entry and separates by indicated symbol 
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the ident ti that column 

# Now we can easily refer to a cell & treatment type, to see the options : 

levels(TreatmentAnnotated) # list the options to perform DE on :

# MYELOID CELL TYPES : 
# 1. "cMono_alpha", "cMono_lambda"           VS    cMono_untreated
# 2. "intMono_alpha" , "intMono_lambda"      VS    intMono_untreated
# 3. "neutrophils_alpha", neutrophils_lambda"VS    neutrophils_untreated 

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

# MYELOID CELL TYPES : 

# 1. Classical monocytes 
M1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cMono_alpha", ident.2 = "cMono_untreated", sep = "_") 
M1AlphaResponse$gene <- rownames(M1AlphaResponse) # puts gene names into a column
cMonoAlphaResponse <- data.frame(Gene = M1AlphaResponse$gene, Log2FoldChange = M1AlphaResponse$avg_log2FC)  
write.csv(cMonoAlphaResponse, "honours/results/DEAnalysis/cMonoAlphaResponse.csv", row.names = FALSE)
# size = 1360 DEGs 

M1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "cMono_lambda", ident.2 = "cMono_untreated", sep = "_") 
M1LambdaResponse$gene <- rownames(M1LambdaResponse) 
cMonoLambdaResponse <- data.frame(Gene = M1LambdaResponse$gene, Log2FoldChange = M1LambdaResponse$avg_log2FC)  
write.csv(cMonoLambdaResponse, "honours/results/DEAnalysis/cMonoLambdaResponse.csv", row.names = FALSE)
# size = 67 DEGs 


# 2. Intermediate monocytes 
M2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "intMono_alpha", ident.2 = "intMono_untreated", sep = "_") 
M2AlphaResponse$gene <- rownames(M2AlphaResponse) 
intMonoAlphaResponse <- data.frame(Gene = M2AlphaResponse$gene, Log2FoldChange = M2AlphaResponse$avg_log2FC)  
write.csv(intMonoAlphaResponse, "honours/results/DEAnalysis/intMonoAlphaResponse.csv", row.names = FALSE)
# size = 1413 DEGs 

M2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "intMono_lambda", ident.2 = "intMono_untreated", sep = "_") 
M2LambdaResponse$gene <- rownames(M2LambdaResponse) 
intMonoLambdaResponse <- data.frame(Gene = M2LambdaResponse$gene, Log2FoldChange = M2LambdaResponse$avg_log2FC)  
write.csv(intMonoLambdaResponse, "honours/results/DEAnalysis/intMonoLambdaResponse.csv", row.names = FALSE)
# size = 286 DEGs 

# 3. Neutrophils 
M3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_alpha", ident.2 = "neutrophils_untreated", sep = "_") 
M3AlphaResponse$gene <- rownames(M3AlphaResponse) 
neuAlphaResponse <- data.frame(Gene = M3AlphaResponse$gene, Log2FoldChange = M3AlphaResponse$avg_log2FC)  
write.csv(neuAlphaResponse, "honours/results/DEAnalysis/neuAlphaResponse.csv", row.names = FALSE)
# size = 1507 DEGs 

M3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "neutrophils_lambda", ident.2 = "neutrophils_untreated", sep = "_") 
M3LambdaResponse$gene <- rownames(M3LambdaResponse) 
neuLambdaResponse <- data.frame(Gene = M3LambdaResponse$gene, Log2FoldChange = M3LambdaResponse$avg_log2FC)  
write.csv(neuLambdaResponse, "honours/results/DEAnalysis/neuLambdaResponse.csv", row.names = FALSE)
# size = 175 DEGs 


# LYMPHOID CELL TYPES : 
# 1. CD4_helper T cells 
L1AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_helper_alpha", ident.2 = "CD4_helper_untreated", sep = "_") 
L1AlphaResponse$gene <- rownames(L1AlphaResponse) 
CD4hAlphaResponse <- data.frame(Gene = L1AlphaResponse$gene, Log2FoldChange = L1AlphaResponse$avg_log2FC)  
write.csv(CD4hAlphaResponse, "honours/results/DEAnalysis/CD4hAlphaResponse.csv", row.names = FALSE)
# size =  762 DEGs 

L1LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_helper_lambda", ident.2 = "CD4_helper_untreated", sep = "_") 
L1LambdaResponse$gene <- rownames(L1LambdaResponse) 
CD4hLambdaResponse <- data.frame(Gene = L1LambdaResponse$gene, Log2FoldChange = L1LambdaResponse$avg_log2FC)  
write.csv(CD4hLambdaResponse, "honours/results/DEAnalysis/CD4hLambdaResponse.csv", row.names = FALSE)
# size =  28 DEGs 

# 2. CD4_naive T cells 
L2AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_naive_alpha", ident.2 = "CD4_naive_untreated", sep = "_") 
L2AlphaResponse$gene <- rownames(L2AlphaResponse) 
CD4nAlphaResponse <- data.frame(Gene = L2AlphaResponse$gene, Log2FoldChange = L2AlphaResponse$avg_log2FC)  
write.csv(CD4nAlphaResponse, "honours/results/DEAnalysis/CD4nAlphaResponse.csv", row.names = FALSE)
# size =  951 DEGs 

L2LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD4_naive_lambda", ident.2 = "CD4_naive_untreated", sep = "_") 
L2LambdaResponse$gene <- rownames(L2LambdaResponse) 
CD4nLambdaResponse <- data.frame(Gene = L2LambdaResponse$gene, Log2FoldChange = L2LambdaResponse$avg_log2FC)  
write.csv(CD4nLambdaResponse, "honours/results/DEAnalysis/CD4nLambdaResponse.csv", row.names = FALSE)
# size =  37 DEGs 

# 3. Tregs_ 
L3AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_alpha", ident.2 = "Tregs_untreated", sep = "_") 
L3AlphaResponse$gene <- rownames(L3AlphaResponse)
TregsAlphaResponse <- data.frame(Gene = L3AlphaResponse$gene, Log2FoldChange = L3AlphaResponse$avg_log2FC)  
write.csv(TregsAlphaResponse, "honours/results/DEAnalysis/TregsAlphaResponse.csv", row.names = FALSE)
# size =  1465 DEGs 

L3LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "Tregs_lambda", ident.2 = "Tregs_untreated", sep = "_") 
L3LambdaResponse$gene <- rownames(L3LambdaResponse) 
TregsLambdaResponse <- data.frame(Gene = L3LambdaResponse$gene, Log2FoldChange = L3LambdaResponse$avg_log2FC)  
write.csv(TregsLambdaResponse, "honours/results/DEAnalysis/TregsLambdaResponse.csv", row.names = FALSE)
# size = 602 DEGs 

# 4. CD8_alpha 
L4AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD8_alpha", ident.2 = "CD8_untreated", sep = "_") 
L4AlphaResponse$gene <- rownames(L4AlphaResponse)
CD8AlphaResponse <- data.frame(Gene = L4AlphaResponse$gene, Log2FoldChange = L4AlphaResponse$avg_log2FC)  
write.csv(CD8AlphaResponse, "honours/results/DEAnalysis/CD8AlphaResponse.csv", row.names = FALSE)
# size =  961 DEGs 

L4LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "CD8_lambda", ident.2 = "CD8_untreated", sep = "_") 
L4LambdaResponse$gene <- rownames(L4LambdaResponse) 
CD8LambdaResponse <- data.frame(Gene = L4LambdaResponse$gene, Log2FoldChange = L4LambdaResponse$avg_log2FC)  
write.csv(CD8LambdaResponse, "honours/results/DEAnalysis/CD8LambdaResponse.csv", row.names = FALSE)
# size = 69 DEGs 

# 5. NK cells 
L5AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_alpha", ident.2 = "NK_untreated", sep = "_") 
L5AlphaResponse$gene <- rownames(L5AlphaResponse)
NKAlphaResponse <- data.frame(Gene = L5AlphaResponse$gene, Log2FoldChange = L5AlphaResponse$avg_log2FC)  
write.csv(NKAlphaResponse, "honours/results/DEAnalysis/NKAlphaResponse.csv", row.names = FALSE)
# size = 2160  DEGs 

L5LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "NK_lambda", ident.2 = "NK_untreated", sep = "_") 
L5LambdaResponse$gene <- rownames(L5LambdaResponse) 
NKLambdaResponse <- data.frame(Gene = L5LambdaResponse$gene, Log2FoldChange = L5LambdaResponse$avg_log2FC)  
write.csv(NKLambdaResponse, "honours/results/DEAnalysis/NKLambdaResponse.csv", row.names = FALSE)
# size = 923 DEGs 

# 6. B CELLS 
L6AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_alpha", ident.2 = "B_untreated", sep = "_") # size = 1287
L6AlphaResponse$gene <- rownames(L6AlphaResponse) # puts gene names into a column
BAlphaResponse <- data.frame(Gene = L6AlphaResponse$gene, Log2FoldChange = L6AlphaResponse$avg_log2FC)  
write.csv(BAlphaResponse, "honours/results/DEAnalysis/BAlphaResponse.csv", row.names = FALSE) # next step requires entry as a csv hence this step 
# size = 1287 DEGs

L6LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_lambda", ident.2 = "B_untreated", sep = "_") # size = 310 
L6LambdaResponse$gene <- rownames(L6LambdaResponse)
BLambdaResponse <- data.frame(Gene = L6LambdaResponse$gene,Log2FoldChange = L6LambdaResponse$avg_log2FC)  
write.csv(BLambdaResponse, "honours/results/DEAnalysis/BLambdaResponse.csv", row.names = FALSE)
# size = 310 DEGs 

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
head(AlphaResponse)
# Create a Mart object for the human genome
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# List available datasets for the Mart 
listDatasets(mart)

# Retrieve gene information for human genes : 
gene_ID <- getBM(attributes = c ("ensembl_gene_id", "external_gene_name", "entrezgene_id"), 
                 mart = mart)
AlphaResponse <- merge(AlphaResponse, gene_ID[,c(2,3)], by.x = "gene", by.y = "external_gene_name")


de <- AlphaResponse$entrezgene_id
edo <- enrichDGN(de) # The enrichDGN function from the DGN package is designed to perform enrichment analysis using the DisGeNET database, which contains gene-disease associations

ego <- enrichGO(gene = de, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05
                )

# Sort enrichment results by p-value
sorted_ego <- ego[order(ego@result$pvalue), ]
topego <- sorted_ego[1:50,]

library(enrichplot)

# Plot the enrichment results
barplot(topego) 
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
