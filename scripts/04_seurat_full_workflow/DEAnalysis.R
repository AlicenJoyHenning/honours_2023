# Differential Expression analysis in R (seurat

###### dependencies & load object #####
library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(patchwork)

treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")

TreatmentAnnotated <- RenameIdents(treatment, 
                                   '0' = 'monocytes1',
                                   '1' = 'monocytes2',
                                   '2' = 'CD4+T_helper_cells',
                                   '3' = 'CD4+_Naive_T_cells',
                                   '4' = 'neutrophils',
                                   '5' = 'CD8+_T_cells',
                                   '6' = 'B_cells',
                                   '7' = 'CD8+_T_cells',
                                   '8' = 'CD8+_T_cells',
                                   '9' = 'NK_cells',
                                   '10' = 'unknown1',
                                   '11' = 'Tregs',
                                   '12' = 'unknown2',
                                   '13' = 'unknown3',
                                   '14' = 'unknown4')

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
