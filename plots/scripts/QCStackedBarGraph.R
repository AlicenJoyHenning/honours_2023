# Plotting 
library(ggplot2)

##### Quality control stacked bar graph  ######

UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz")
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)

AlphaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/matrix.mtx.gz","honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)

LambdaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)


dim(alphamg)
dim(lambdamg)
dim(untreatedmg)

alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

alphalg <- subset(alpha, subset = nFeature_RNA < 200)
alphahg <- subset(alpha, subset = nFeature_RNA > 2500)
alphamg <- subset(alpha, subset = percent.mt > 10)
alphamg <- subset(alpha, subset = nFeature_RNA >= 200 & nFeature_RNA <= 2500)
alphamg

lambdalg <- subset(lambda, subset = nFeature_RNA < 200)
lambdahg <- subset(lambda, subset = nFeature_RNA > 2500)
lambdamg <- subset(lambda, subset = percent.mt > 10)
lambdamg

untreatedlg <- subset(untreated, subset = nFeature_RNA > 200)
untreatedhg <- subset(untreated, subset = nFeature_RNA < 2500)
untreatedmg <- subset(untreated, subset = percent.mt < 10)






QCMetric <- rep(c("LOW nF", "Percent MT", "HIGH nF"), 3)
cells <- c(0, 6169, 331, 5838, 1079, 5090)
grouped <- data.frame(QCMetric, cells, colours)
colours <- c("LOW nF" = "#6ab5ba", "Percent Mt"="#6ab5ba", "HIGH nF"="#6ab5ba")
labels <- c("0%", "5%", "18%")


vertical <- 
  ggplot(grouped,
  aes(fill=colours, y=cells, x=QCMetric)) + 
  geom_bar(color = "black", position="stack", stat="identity") + 
  geom_text(inherit.aes = TRUE, aes(label= labels), vjust=0) +
 geom_label(label = labels, aes(fill = QCMetric), colour = "black", fontface = "bold") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = colours)+ 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12, face = "bold"),   
    axis.text.y = element_text(size = 12),                 
    axis.title = element_text(size = 14, face = "bold"),   
    legend.title = element_blank(),                        
    legend.text = element_text(size = 12)                
  )





horizontal <-   
  ggplot(
  grouped,
  aes(fill=QCMetric, y=treatments, x=cells)) + 
  geom_bar(color = "black", position="stack", stat="identity") + 
  geom_text(inherit.aes = TRUE, aes(label= labels), vjust=0) +
  geom_label(label = labels, aes(fill = QCMetric), colour = "black", fontface = "bold") +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = colours)+ 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12, face = "bold"),   # Adjust x-axis text size and style
    axis.text.y = element_text(size = 12),                 # Adjust y-axis text size
    axis.title = element_text(size = 14, face = "bold"),   # Adjust axis title size and style
    legend.title = element_text(size = 14, face = "bold"),                        # Remove legend title
    legend.text = element_text(size = 12)                  # Adjust legend text size
  )



colours <- c("LOW nF" = "#6ab5ba", "Percent Mt" = "#6ab5ba", "HIGH nF" = "#6ab5ba")
labels <- c("0%", "5%", "18%")

grouped <- data.frame(QCMetric = c("LOW nF", "Percent Mt", "HIGH nF"),
                      cells = c(0, 6169, 331, 5838, 1079, 5090))

vertical <- 
  ggplot(grouped,
         aes(fill = QCMetric, y = cells, x = QCMetric)) + 
  geom_bar(color = "black", position = "stack", stat = "identity") + 
  geom_text(inherit.aes = TRUE, aes(label = labels), vjust = 0) +
  geom_label(aes(label = labels), colour = "black", fontface = "bold", position = "stack", vjust = 0) +
  theme(panel.background = element_rect(fill = "white")) +
  scale_fill_manual(values = colours) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = 12, face = "bold"),   
    axis.text.y = element_text(size = 12),                 
    axis.title = element_text(size = 14, face = "bold"),   
    legend.title = element_blank(),                        
    legend.text = element_text(size = 12)
  )

print(vertical)
  


labels <- c("0%", "5%", "18%")

grouped <- data.frame(QCMetric = c("LOW nF", "Percent Mt", "HIGH nF"),
                      cells = c(6169, 5838, 5090),
                      stack = c(0, 331, 1079))
vertical <- 
  ggplot(grouped,
         aes(fill = QCMetric, y = cells, x = QCMetric)) +
  geom_bar(color = "black", position = "stack", stat = "identity") + 
  geom_bar(aes(y = stack), color = "black", fill = "#6ab5ba", position = "stack", stat = "identity") +
  geom_text(aes(label = labels), position = position_stack(vjust = 0.5)) +
  geom_label(aes(label = labels), colour = "black", fontface = "bold", position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),   
    axis.text.y = element_text(size = 12),                 
    axis.title = element_text(size = 14, face = "bold"),   
    legend.title = element_blank(),                        
    legend.text = element_text(size = 12)
  ) +
  labs(title = "Stacked Bar Graph", x = "QCMetric", y = "Number of Cells")

print(vertical)









