# Plotting the filtering of cells to work with based on QC metrics 

##### Loading dependencies  ######
# Plotting 
library(ggplot2)
library(scales)

##### Quality control stacked bar graph : data prep ######

UntreatedMatrix <- ReadMtx("honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/untreatedDarisiaIndex/seurat_matrix/features.tsv.gz")
untreated <- CreateSeuratObject(UntreatedMatrix, project="untreated", min.cells=3, min.features=200)

AlphaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/matrix.mtx.gz","honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnalphaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
alpha <- CreateSeuratObject(AlphaMatrix, project="alpha", min.cells=3, min.features=200)

LambdaMatrix <- ReadMtx("honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/matrix.mtx.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/barcodes.tsv.gz", "honours/work/DarisiaIndex/ifnlambdaDarisiaIndex/seurat_matrix/AdjustedFeatures.tsv.gz", skip.feature = 1)
lambda <- CreateSeuratObject(LambdaMatrix, project="lambda", min.cells=3, min.features=200)


dim(lambdamg)
dim(untreatedmg)

alpha[["percent.mt"]] <- PercentageFeatureSet(alpha, pattern = "^MT-")
lambda[["percent.mt"]] <- PercentageFeatureSet(lambda, pattern = "^MT-")
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

alphalg <- subset(alpha, subset = nFeature_RNA <= 200 & percent.mt < 10)
alphahg <- subset(alpha, subset = nFeature_RNA >= 2500 & percent.mt < 10)
alphamg <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt >= 10)
alphahm <- subset(alpha, subset = nFeature_RNA >= 2500 & percent.mt >= 10)
alphalm <- subset(alpha, subset = nFeature_RNA <= 200 & percent.mt >= 10)
alphaf <- subset(alpha, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
dim(alphahm)

lambdalg <- subset(lambda, subset = nFeature_RNA <= 200 & percent.mt < 10)
lambdahg <- subset(lambda, subset = nFeature_RNA >= 2500 & percent.mt < 10)
lambdamg <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt >= 10)
lambdahm <- subset(lambda, subset = nFeature_RNA >= 2500 & percent.mt >= 10)
lambdalm <- subset(lambda, subset = nFeature_RNA <= 200 & percent.mt >= 10)
lambdaf <- subset(lambda, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
dim(lambdaf)


untreatedlg <- subset(untreated, subset = nFeature_RNA <= 200 & percent.mt < 10)
untreatedhg <- subset(untreated, subset = nFeature_RNA >= 2500 & percent.mt < 10)
untreatedmg <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt >= 10)
untreatdhm <- subset(untreated, subset = nFeature_RNA >= 2500 & percent.mt >= 10)
untreatedlm <- subset(untreated, subset = nFeature_RNA <= 200 & percent.mt >= 10)
untreatedf <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
dim(untreatedlm)


##### Plotting Altogether #####

treatments <- c(rep("alpha", 4), rep("lambda", 4), rep("untreated", 4))
QCMetric <- rep(c("LOW nF", "Percent MT", "HIGH nF","Accepted"), 3)
cells <- c(0, 330, 1072, 4767, 0, 476, 947, 4511,0, 399, 1430, 4656)
grouped <- data.frame(treatments, QCMetric, cells)
labels <- c("0%", "5%", "18%", "77%","0%", "8%", "16%", "76%","0%", "6%", "22%", "72%")
colours <- c("LOW nF" = "#289ce5", "Percent Mt"="#06b96b", "HIGH nF" ="#a4dedd", "Accepted" = "#addfbb")

vertical <- 
  ggplot(
    grouped,
    aes(fill=QCMetric, y=cells, x=treatments)) + 
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
    legend.title = element_blank(),                        # Remove legend title
    legend.text = element_text(size = 12)                  # Adjust legend text size
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

##### Plotting sets individually ####

palette.d <- c("#7c8c94", # high
               "#6ab5ba", # high m
               "#c9cacd", # kept
              #"#a1a1a1", # low 
               "#696969", # low m 
               "#9dcfd3") # m

# ALPHA : 
# Create data frame 
alphapie <- data.frame(
  Metric=c("kept","low count\n& mt percent", "high count","mt percent", "high count\n& mt percent"),
  count=c(3160, 3, 1084, 1752, 210)
)

# Compute percentages
alphapie$fraction <- alphapie$count / sum(alphapie$count)

# Compute the cumulative percentages (top of each rectangle)
alphapie$ymax <- cumsum(alphapie$fraction)

# Compute the bottom of each rectangle
alphapie$ymin <- c(0, head(alphapie$ymax, n=-1))

# Compute label position
alphapie$labelPosition <- (alphapie$ymax + alphapie$ymin) / 2

# Compute percentages for labels 
alphapie$label <- percent(alphapie$fraction, accuracy = 1, suffix = "")  # Use accuracy = 1 to remove decimal points # suffix to remove percent sign 

# Make the plot
alpha <- ggplot(alphapie, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Metric)) +
    geom_rect() + 
# geom_rect(color = "black", size = 0.65) +
# geom_label( x=3.5, aes(y=labelPosition, label=label), size=8) +
  geom_text(aes(x = 3.75, y = labelPosition, label = label,  fontface = "bold"), size = 4) +  # Add labels with percentages
  scale_fill_manual(values = palette.d)+ 
  coord_polar(theta="y") +
  xlim(c(2.7, 4)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"))


# LAMBDA : 
# Create data frame 
lambdapie <- data.frame(
  category=c("kept","high count","low count\n& mt percent", "mt percent", "high count\n& mt percent"),
  count=c(3020, 1253, 3, 1975, 455)
)

# Compute percentages
lambdapie$fraction <- lambdapie$count / sum(lambdapie$count)

# Compute the cumulative percentages (top of each rectangle)
lambdapie$ymax <- cumsum(lambdapie$fraction)

# Compute the bottom of each rectangle
lambdapie$ymin <- c(0, head(lambdapie$ymax, n=-1))

# Compute label position
lambdapie$labelPosition <- (lambdapie$ymax + lambdapie$ymin) / 2

# Compute percentages for labels 
lambdapie$label <- percent(lambdapie$fraction, accuracy = 1, suffix = "")  # Use accuracy = 1 to remove decimal points

# Make the plot
lambda <- ggplot(lambdapie, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect() + 
  # geom_rect(color = "black", size = 0.65) +
  # geom_label( x=3.5, aes(y=labelPosition, label=label), size=8) +
  geom_text(aes(x = 3.75, y = labelPosition, label = label,  fontface = "bold"), size = 4) +
  scale_fill_manual(values = palette.d)+ 
  coord_polar(theta="y") +
  xlim(c(2.7, 4)) +
  theme_void() +
  theme(legend.position = "none")


# UNTREATED : 
# Create data frame 
untreatedpie <- data.frame(
  QCmetric =c("kept", "high count",  "low count\n& mt percent", "mt percent", "high count\n& mt percent"),
  count=c(2870, 782, 5, 2243, 336)
)

# Compute percentages
untreatedpie$fraction <- untreatedpie$count / sum(untreatedpie$count)

# Compute the cumulative percentages (top of each rectangle)
untreatedpie$ymax <- cumsum(untreatedpie$fraction)

# Compute the bottom of each rectangle
untreatedpie$ymin <- c(0, head(untreatedpie$ymax, n=-1))

# Compute label position
untreatedpie$labelPosition <- (untreatedpie$ymax + untreatedpie$ymin) / 2

# Compute percentages for labels 
untreatedpie$label <- percent(untreatedpie$fraction, accuracy = 1, suffix = "")  # Use accuracy = 1 to remove decimal points

# Make the plot
untreated <- ggplot(untreatedpie, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=QCmetric)) +
  geom_rect() + 
  # geom_rect(color = "black", size = 0.65) +
  # geom_label( x=3.5, aes(y=labelPosition, label=label), size=8) +
  geom_text(aes(x = 3.75, y = labelPosition, label = label,  fontface = "bold"), size = 4) +
  scale_fill_manual(values = palette.d)+ 
  coord_polar(theta="y") +
  xlim(c(2.7, 4)) +
  theme_void() +
  theme(legend.position = "none")



# all 

untreated | lambda | alpha


