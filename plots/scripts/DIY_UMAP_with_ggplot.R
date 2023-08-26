# This script is for making (customizing a UMAP plot using ggplot :
library(patchwork)
library(gridExtra)
##### Overall Plot ####

treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")

# Extract UMAP coordinatesb and cluster information
treatment.umap.coords <- as.data.frame(treatment@reductions$umap@cell.embeddings)
clusters <- treatment$seurat_clusters

# Create a dataframe for ggplot
treatment.df <- data.frame(
  x = treatment.umap.coords$UMAP_1,
  y = treatment.umap.coords$UMAP_2,
  seurat_clusters = factor(clusters)
)

# Define color palette
palette.a <- RColorBrewer::brewer.pal(12, "Paired")
palette.b <- c("#ee5e17", #0
               "#d72554", #1
               "#6ab5ba", #2
               "#2e8f95", #3
               "#900c3e", #4
               "#8caf2e", #5
               "#b4d56c", #6
               "#297b57", #7 
               "#00945a", #8
               "#265221", #9
               "#FFCB3E", #10
               "#00a68e", #11
               "#5c040c", #12
               "#f7bc6e", #13
               "#6ab5ba" #14
              )

# Create the ggplot plot
total <- 
  ggplot(treatment.df, aes(x, y, colour = seurat_clusters)) +
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



##### Colour each cluster individually #####

# Select the cluster you want to color
selected_cluster <- 5  # Change this to the desired cluster number



##### NK cell type : cluster 9  #####
NKCluster <- 10
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[NKCluster] <- palette.b[NKCluster]

NKLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )


features = c("CD96", "APOBEC3G", "CTSW","NKG7", "GNLY", "GZMA", "FCGR3A")

NKFP <- 
  FeaturePlot(object = treatment, 
            features = "NKG7",
            cols = c("lightgrey", "black"),
            label = TRUE,
            pt.size = 1.5, 
            blend = FALSE, 
            interactive = FALSE) + theme(
              panel.background = element_rect(fill = "darkgrey")) 

cluster_boundaries <- c(5, 7)

NKDP <-
  DotPlot(object = treatment, 
        features = features,
        cols = c("grey", "#265221")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Adjust x-axis labels angle
  geom_hline(yintercept = c(9, 11), linetype = "dotted", color = "black")

NK <- (NKLook / NKFP) | NKDP # uses patchwork since ggplot2 was not used to make the feature plot and Lord knows I'm not going to use it 

##### B cell type : cluster 6 #####
BCluster <- 7
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[BCluster] <- palette.b[BCluster]

BLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )
##### Reguatory T cells #####
TregsCluster <- 12
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[TregsCluster] <- palette.b[TregsCluster]

TregsLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )


features = c("FOXP3", "TIGIT", "TRAC", "IL32")

TregFP <- 
  FeaturePlot(object = treatment, 
              features = "FOXP3",
              cols = c("lightgrey", "black"),
              label = TRUE,
              pt.size = 1.5, 
              blend = FALSE, 
              interactive = FALSE) + theme(
                panel.background = element_rect(fill = "darkgrey")) 

TregDP <-
  DotPlot(object = treatment, 
          features = features,
          cols = c("grey", "#00a68e")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(11, 13), linetype = "dotted", color = "black")

Treg <- (TregsLook / TregFP) | TregDP
Treg

##### CD8T cells #####
CD8Cluster <- c(6, 8, 9)
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[CD8Cluster] <- palette.b[CD8Cluster]

CD8Look <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### Naive T @TregsCluster <- 12
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[TregsCluster] <- palette.b[TregsCluster]

TregsLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )





##### Naive T #####
NaiveTCluster <- 4
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[NaiveTCluster] <- palette.b[NaiveTCluster]

NaiveTLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### CD4 T  ####
CD4Cluster <- 3
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[CD4Cluster] <- palette.b[CD4Cluster]

CD4Look <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### Neutrophils #####
NCluster <- 5
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[NCluster] <- palette.b[NCluster]

NLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### Monocytes ####
MCluster <- c(1,2)
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[MCluster] <- palette.b[MCluster]

MLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### Altogether #####
# combined_layout <- (
#   total |
#     (
#       (NKLook | CD8Look | NaiveLook | BLook) / 
#         (TregsLook | CD4Look | NLook | MLook)
#     )
# )

combined_layout <- (
    
      (NKLook | CD8Look | NaiveLook | BLook) / 
        (TregsLook | CD4Look | NLook | MLook)
    
)

combined_plot <- combined_layout + 
  plot_layout(guides = "collect")
