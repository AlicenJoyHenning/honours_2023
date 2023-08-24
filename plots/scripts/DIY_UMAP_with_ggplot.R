# This script is for making (customizing a UMAP plot using ggplot :

# Alternatively, using ggplot :   
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
palette.b <- c(palette.a,"#2E8B57","#7FCDBB","#41B6C4","#1D91C0")

# Create the ggplot plot
ggplot(treatment.df, aes(x, y, colour = seurat_clusters)) +
  geom_point(size = 1) +
  scale_colour_manual(values = palette.b) +
  labs(#title = "IFN alpha",
    x = "UMAP 1",  # Rename x-axis label
    y = "UMAP 2",
    color = "")  + 
  theme_minimal() + 
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
      fill = palette.a, 
      size = 3.5),  # Color the squares with the same palette
    key_height = unit(1, "npc"),  # Spread the legend dots across the vertical length
    key_width = unit(4, "cm"),   # Adjust the width of the legend blocks
    title.theme = element_text(hjust = 0.5),  # Center the legend title
    label.position = "right",
    label.hjust = 1
  ))
