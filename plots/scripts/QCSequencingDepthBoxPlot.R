# Sequencing Depth Boxplot 

library(Seurat) # using Seurat object 
library(Matrix)
library(ggplot2)
library(dplyr)

setwd("")
getwd()


##### Load datasets: alpha, lambda, and untreated using ReadMtx function :  ####
Amatrix <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/matrix.mtx"
Abarcodes <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/unfiltered_counts.barcodes.txt"
Afeatures <- "D:/AlicenHonours/work/1109/alpha/Count_analysis/unfiltered_counts/unfiltered_counts.genes.txt"
AlphaMatrix <- ReadMtx(Amatrix, Abarcodes, Afeatures)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=0)

Lmatrix <- "D:/AlicenHonours/work/1109/lambda/Count_analysis/unfiltered_counts/unfiltered_counts.mtx"
Lbarcodes <- "D:/AlicenHonours/work/1109/lambda/Count_analysis/unfiltered_counts/unfiltered_counts.barcodes.txt"
Lfeatures <- "D:/AlicenHonours/work/1109/lambda/Count_analysis/unfiltered_counts/unfiltered_counts.genes.txt"
LambdaMatrix <- ReadMtx(Lmatrix, Lbarcodes, Lfeatures)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features = 0, feature.colum)

Umatrix <- "honours/work/1109/untreated/matrix.mtx.gz"
Ubarcodes <- "honours/work/1109/untreated/barcodes.tsv.gz"
Ufeatures <- "honours/work/1109/untreated/features.tsv.gz"
UntreatedMatrix <- ReadMtx(Umatrix, Ubarcodes, Ufeatures)
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features = 0)



##### Boxplot ####
# loading the datasets into data frames for plotting: 
alpha <- as.data.frame(alpha@meta.data, 
                       orig.ident = alpha@meta.data$orig.ident,
                       nCount_RNA = alpha@meta.data$nCount_RNA,
                       nFeature_RNA = alpha@meta.data$nFeature_RNA)
lambda <- as.data.frame(lambda@meta.data, 
                       orig.ident = lambda@meta.data$orig.ident,
                       nCount_RNA = lambda@meta.data$nCount_RNA,
                       nFeature_RNA = lambda@meta.data$nFeature_RNA)
untreated <- as.data.frame(untreated@meta.data, 
                       orig.ident = untreated@meta.data$orig.ident,
                       nCount_RNA = untreated@meta.data$nCount_RNA,
                       nFeature_RNA = untreated@meta.data$nFeature_RNA)
# combining the datasets :
combined <- rbind(alpha, lambda, untreated)

combinedalt <- subset(combined, nCount_RNA <= 15000)
combined2 <- subset(combined, nCount_RNA >= 15000)
combined <- subset(combined, nCount_RNA <= 50000)

# Labeling outliers in metadata column 
combined$outlier <- with(combined, ifelse(nCount_RNA >= 15000 & nCount_RNA <= 50000, TRUE, FALSE))



# lower portion of box plot 
# box <- ggplot(combinedalt, 
#                aes(x = orig.ident, y = nCount_RNA)) + 
#   geom_boxplot() + 
#   theme_classic()
# 
# 
# 
# jitter <- ggplot(combined2, 
#                aes(x = orig.ident, y = nCount_RNA)) + 
#   geom_jitter(alpha = 0.5, width = 0.2, color = "grey", shape = 16) +  # Add jitter for better visibility of outliers
#   theme_classic()


combined_plot <- ggplot(combined, aes(x = orig.ident, y = nCount_RNA, fill = outlier)) +
  geom_boxplot(data = subset(combined, outlier == FALSE), fill = "grey") +
  geom_jitter(data = subset(combined, outlier == TRUE), width = 0.18, color = "grey", shape = 16) +
  theme_classic() +
  labs(x = "", y = "Total UMI count per barcode") +
  guides(color = FALSE) +  # Remove the legend
  scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "red")) + # Customize box colors
  scale_x_discrete(labels = c("u", "λ", "α")) +
  theme(axis.text.x = element_text(size = 14, face =  "bold"),  
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 14, face =  "bold"),
        panel.border = element_rect(color = "white", fill = NA, linewidth  = 1),
        axis.line.x = element_blank(),  # Remove x-axis
        axis.line.y = element_blank())  # Remove y-axis  # Add black border around the plot area)  


# 
# box <- ggplot(combined, aes(x = orig.ident, y = nCount_RNA)) +
#   geom_boxplot() +
#   geom_jitter(alpha = 0.5, width = 0.2, color = "grey", shape = 16) +  # Add jitter for better visibility of outliers
#   theme_classic()




