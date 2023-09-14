# Cluster QC Boxplots 
library(ggplot2)
library(patchwork)

##### Data set generation ####

treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")
# Create a dataframe for ggplot

clusterQCdf <- data.frame(
  count = treatment$nCount_RNA,
  feature = treatment$nFeature_RNA,
  mito = treatment$percent.mt,
  clusters = treatment$seurat_clusters
)

doublet <- data.frame(
  alpha = alpha@meta.data$percent.xist
)

clustercount <- 
  ggplot(alpha, aes(x=doublet, y=alpha, fill=doublet)) + 
  geom_boxplot(color="black", fill="#6ab5ba",alpha = 0.7) + #alpha=0.3
  theme(legend.position="none") 

##### Generate the box plot #####

clustercount <- 
  ggplot(clusterQCdf, aes(x=clusters, y=count, fill=clusters)) + 
  geom_boxplot(color="black", fill="#6ab5ba",alpha = 0.7) + #alpha=0.3
  theme(legend.position="none") 


clusterfeature <- 
  ggplot(clusterQCdf, aes(x=clusters, y=feature, fill=clusters)) + 
  geom_boxplot(color="black", fill="#6ab5ba",alpha = 0.7) + #alpha=0.3
  theme(legend.position="none") 

##### Pie chart for cells  #####

clusterQCdf$cells <- rownames(clusterQCdf)

# Count the number of identical entries in the "Name" column
cells <- table(clusterQCdf$clusters)
clustercount <- data.frame(cells = cells)
names(clustercount) <- c("cluster", "freq")
clustercount$cluster@0

# counting the cells that are in clusters being used 
IdentifiedCells <- clustercount[clustercount$cluster %in% c(0:9, 11),]  
IdentifiedCellsSum <- sum(IdentifiedCells$freq)

# counting the cells that are in clusters being NOT being used 
UnidentifiedCells <- clustercount[clustercount$cluster %in% c(10, 12, 13, 14),]  
UnidentifiedCellsSum <- sum(UnidentifiedCells$freq)

# Percentages 
IdentifiedPer <- ((IdentifiedCellsSum)/(IdentifiedCellsSum + UnidentifiedCellsSum))* 100
UnidentifiedPer <- ((UnidentifiedCellsSum)/(IdentifiedCellsSum + UnidentifiedCellsSum))* 100

# creating data frame for pie chart 
pieData <- c(IdentifiedCellsSum, UnidentifiedCellsSum)

labels <- c("identified\ncells\n95.98%", "unidentified\ncells\n4.02%")
pie <- pie(pieData, labels = labels, border = "black" ,col = c("#6ab4b9", "grey"), radius = 2, lwd = 4)

###### Create the pie chart using ggplot2 #####
# Sample data
IdentifiedCellsSumP <- 95.98
UnidentifiedCellsSumP <- 4.02

# Create a data frame for the pie chart
pie_df <- data.frame(
  label = c("identified\ncells", "unidentified\ncells"),
  value = c(IdentifiedCellsSumP, UnidentifiedCellsSumP)
)

# Create the pie chart using ggplot2
pie_chart <- ggplot(pie_df, aes(x = "", y = value, fill = label)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = c("#6ab5ba", "grey")) +
  labs(fill = NULL) +
  theme(legend.position = "none") +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank())

# Print the pie chart
print(pie_chart)

##### Overall #####
clustercount
