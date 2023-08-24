# This script is for making (customizing a UMAP plot using ggplot :

# Alternatively, using ggplot :   

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
palette.b <- c("#FB836F", #0
               "#d72554", #1
               "#6ab5ba", #2
               "#2e8f95", #3
               "#7E549F", #4
               "#8caf2e", #5
               "#69a923", #6
               "#297b57", #7 
               "#00945a", #8
               "#265221", #9
               "#FFCB3E", #10
               "#00a68e", #11
               "#5c040c", #12
               "#ef931b", #13
               "#6ab5ba" #14
              )

# Create the ggplot plot
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
selected_cluster <- 7  # Change this to the desired cluster number

# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[selected_cluster] <- palette.b[selected_cluster]

ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = ""
  ) + 
  theme_classic() + 
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
