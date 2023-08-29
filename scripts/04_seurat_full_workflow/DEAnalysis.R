# Differential Expression analysis in R (seurat

###### dependencies & load object #####
library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(patchwork)

treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")

TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'class_monocytes',
                                   '1' = 'inter_monocytes',
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

##### View the data #####
cluster_labels <- Idents(TreatmentAnnotated)

# Convert cluster_labels to a dataframe
cluster_labels_df <- data.frame(cell_id = rownames(TreatmentAnnotated@meta.data),
                                cluster_label = cluster_labels)



# Merge cluster_labels_df with TreatmentAnnotated.df based on cell_id
TreatmentAnnotated.df <- merge(TreatmentAnnotated.df, cluster_labels_df,
                               by = "", all.x = TRUE)

###
Idents(TreatmentAnnotated)

TreatmentAnnotated.umap.coords <- as.data.frame(TreatmentAnnotated@reductions$umap@cell.embeddings)
clusters <- TreatmentAnnotated$seurat_clusters

# Create a dataframe for ggplot
TreatmentAnnotated.df <- data.frame(
  x = TreatmentAnnotated.umap.coords$UMAP_1,
  y = TreatmentAnnotated.umap.coords$UMAP_2,
  seurat_clusters = factor(clusters)
)

ggplot(TreatmentAnnotated.df, aes(x, y, colour = seurat_clusters)) +
  geom_point(size = 1) +
  scale_colour_manual(values = palette.b) +
  labs(#title = "IFN alpha",
    x = "UMAP 1",  # Rename x-axis label
    y = "UMAP 2",
    color = "")  + 
  theme_classic() + 
  theme(#panel.background = element_rect(fill = "lightgrey"),  # Set background color
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),  # Increase axis label size
    axis.title = element_text(size = 14), 
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    #legend.background = element_rect(color = "black", fill = "white"),
    legend.position = "right", 
    legend.title = element_text(size =  14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))) + 
  guides(color = guide_legend(
    override.aes = list(
      #hape = rep(22, length(palette.a)),  # Use squares (blocks)
      fill = palette.b, 
      size = 3.5),  # Color the squares with the same palette
    key_height = unit(1, "npc"),  # Spread the legend dots across the vertical length
    key_width = unit(4, "cm"),   # Adjust the width of the legend blocks
    title.theme = element_text(hjust = 0.5),  # Center the legend title
    label.position = "right",
    label.hjust = 1
  ))


##### DE ######
levels(TreatmentAnnotated) # list the options to perform DE on 
# [1] "monocytes1"          "monocytes2"          "CD4+ T helper cells" "Naive CD4+ T cells"  "neutrophils"         "CD8+ T cells"       
# [7] "B cells"             "NK cells"            "unknown1"            "Tregs"               "unknown2"            "unknown3"           
# [13] "unknown4" 

# creating metadata column to hold cell type and stimulation information
TreatmentAnnotated$celltype.stim <- paste(Idents(TreatmentAnnotated), TreatmentAnnotated$stim, sep = "_")
TreatmentAnnotated$celltype <- Idents(TreatmentAnnotated) # restores the cell type column 
Idents(TreatmentAnnotated) <- "celltype.stim" # switch the ident ti that column 

# Use FindMarkers() to find the genes that are different between stimulated and untreated cell types
AlphaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_cells_alpha", ident.2 = "B_cells_untreated", sep = "_")
# size = 1287
AlphaResponse$gene <- rownames(AlphaResponse) # put gene names into a column
BCellAlphaResponse <- data.frame(
  Gene = AlphaResponse$gene,
  Log2FoldChange = AlphaResponse$avg_log2FC
)
write.csv(BCellAlphaResponse, "honours/results/DEAnalysis/BCellAlphaResponse.csv", row.names = FALSE)



# Same for those different upon IFN lambda treatment 
LambdaResponse <- FindMarkers(TreatmentAnnotated, ident.1 = "B_cells_lambda", ident.2 = "B_cells_untreated", sep = "_")
# size = 310 
LambdaResponse$gene <- rownames(LambdaResponse) # put gene names into a column 




# To find common genes : 
B_cells_common_genes <- intersect(AlphaResponse$gene, LambdaResponse$gene) # identify common genes
# Extract common genes from each dataframe 
B_cells_common_alpha <- AlphaResponse[AlphaResponse$gene %in% B_cells_common_genes, ]
B_cells_common_lambda <- LambdaResponse[LambdaResponse$gene %in% B_cells_common_genes, ]

# merge together the 2 subsetted dataframes to create 1 containing only common genes : 
B_cells_common_dataframe <- merge(B_cells_common_alpha, B_cells_common_lambda, by = 'gene', all = TRUE)
# size = 223

B_cells_alpha_genes <- AlphaResponse$gene[!(AlphaResponse$gene %in% LambdaResponse$gene)] # - (all alpha in lambda) = all alpha not in lambda
B_cells_alpha <- AlphaResponse[AlphaResponse$gene %in% B_cells_alpha_genes, ]
# size = 1064 genes unique to alpha 

B_cells_lambda_genes <- LambdaResponse$gene[!(LambdaResponse$gene %in% AlphaResponse$gene)]
B_cells_lambda <- LambdaResponse[LambdaResponse$gene %in% B_cells_lambda_genes, ] 
# size = 87 genes unique to lambda

B_cells <- list('B_cells_Alpha_Response' = AlphaResponse, 'B_cells_Lambda_Response' = LambdaResponse,'B_cells_common_response' = B_cells_common_dataframe, 'B_cells_alpha_specific' = B_cells_alpha,'B_cells_lambda_specific' = B_cells_lambda)
openxlsx::write.xlsx(B_cells, file = "honours/results/DEAnalysis/B_cells.xlsx")

##### Visualisations ##### 

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
BiocManager::install("clusterProfiler", version = "3.14")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(biomaRt)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
# to identify which biological processes or pathways are enriched within the differentially expressed genes

?clusterProfiler # not helpful 
ls("package:clusterProfiler")
# use the 

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
