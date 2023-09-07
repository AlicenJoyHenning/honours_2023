
# CUSTOMIZING UMAP PLOT USING GGPLOT 

##### [0] Dependencies #####
library(patchwork)
library(gridExtra)
library(ggplot2)
library(Seurat)
library(SeuratObject)

##### [1] Overall UMAP plot #####
# Load the data : 
getwd()
setwd("C:/Users/alice")
treatment <- readRDS("honours/results/IntegratedMarkers/treatment.rds")

# Extract UMAP coordinates and cluster information : 
treatment.umap.coords <- as.data.frame(treatment@reductions$umap@cell.embeddings)
clusters <- treatment$seurat_clusters


# Create a data frame for ggplot : 
treatment.df <- data.frame(
  x = treatment.umap.coords$UMAP_1,
  y = treatment.umap.coords$UMAP_2,
  seurat_clusters = factor(clusters) # cluster numbers 
)

# Define color palette : 
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

# Create the ggplot plot : 
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


##### [2] Treatment-specific UMAP plots #####

clusters.new <- treatment$stim

# Create a dataframe for ggplot : 
treatment.df.new <- data.frame(
  x = treatment.umap.coords$UMAP_1,
  y = treatment.umap.coords$UMAP_2,
  clusters = factor(clusters.new)
)

colours <- c("#c35cad","#6ab5ba","#d3d3d3")

# Define custom colors based on the 'stim' column : 
alpha_cluster_color <- ifelse(treatment.df.new$clusters == "alpha", colours[1], "grey")

alpha.plot <- 
  ggplot(treatment.df.new, aes(x, y, colour = clusters)) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = c("#c35cad", "#d3d3d3", "#d3d3d3")) +
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
    legend.position = "none", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

print(alpha.plot)

# Define custom colors based on the 'stim' column: 
lambda_cluster_color <- ifelse(treatment.df.new$clusters == "lambda", colours[2], "grey")

lambda.plot <- 
  ggplot(treatment.df.new, aes(x, y, colour = clusters)) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = c("#d3d3d3", "#6ab5ba", "#d3d3d3")) +
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
    legend.position = "none", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

print(lambda.plot)

# Define custom colors based on the 'stim' column:
untreated_cluster_color <- ifelse(treatment.df.new$clusters == "untreated", colours[3], "grey")

untreated.plot <- 
  ggplot(treatment.df.new, aes(x, y, colour = clusters)) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = c("white","white", "black")) +
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
    legend.position = "none", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

print(untreated.plot)


##### [3.1] NK cell type   #####
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
            features = "FCGR1A",
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

##### [3.2] B cell type  #####
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


features = c("CD79A", "BCL11A", "BANK1", "BCL2", "CD40", "CD19", "BLNK")

BFP <- 
  FeaturePlot(object = treatment, 
              features = "CD79A",
              cols = c("lightgrey", "black"),
              label = TRUE,
              pt.size = 1.5, 
              blend = FALSE, 
              interactive = FALSE) + theme(
                panel.background = element_rect(fill = "darkgrey")) 

BDP <-
  DotPlot(object = treatment, 
          features = features,
          cols = c("grey", "#c35cad")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(6, 8), linetype = "dotted", color = "black")

B <- (BLook / BFP) | BDP
B

##### [3.3] Regulatory T cells #####
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


features = c("TBC1D4","TIGIT", "FOXP3", "IKZF2", "RTKN2", "IL2RA", "PDP1")

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(11, 13), linetype = "dotted", color = "black")+
  geom_hline(yintercept = c(0), linetype = "solid", color = "black")
  

Treg <- (TregsLook / TregFP) | TregDP
Treg

##### [3.4] CD8T cells #####
CD8Cluster <- c(6, 8, 9)
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[CD8Cluster] <- palette.b[CD8Cluster]

CD8TLook <- 
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


features = c("CD8B", "CD8A", "CCR6", "CTSW", "CXCR6", "C1orf21", "GZMK", "GZMA","GZMB", "GZMH", "NKG7", "SAMD3", "FCGR3A", "FGFBP2", "KLRF1", "KLRD1", "GNLY")
features =c("CD8B", "CXCR6", "GZMK", "CD8A", "CTSW", "NKG7", "SAMD3", "KLRD1", "GNLY", "GZMA", "C1orf21", "GZMH", "FGFBP2", "GZMB", "KLRF1") 

CD8TFP <- 
  FeaturePlot(object = treatment, 
              features = "CD8B",
              cols = c("lightgrey", "black"),
              label = TRUE,
              pt.size = 1.5, 
              blend = FALSE, 
              interactive = FALSE) + theme(
                panel.background = element_rect(fill = "darkgrey")) 

CD8TDP <-
  DotPlot(object = treatment, 
          features = features,
          cols = c("grey", "#8caf2e")) +  
  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(5.75, 6.25, 7.5, 10.5), linetype = "dotted", color = "black")


CD8T <- (CD8TLook /CD8TFP) | CD8TDP | (NKLook / NKFP) 
CD8T


##### [3.5] CD4 :  Helper + Naive ####
HelpTCluster <- c(3, 4)
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[HelpTCluster] <- palette.b[HelpTCluster]

HelpTLook <- 
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


features = c("FHIT","CCR4","TSHZ2","CD4", "MAL" ,"TCF7", "KLF2",  "IL7R")

HelpTFP <- 
  FeaturePlot(object = treatment, 
              features = "CD2",
              cols = c("lightgrey", "black"),
              label = TRUE,
              pt.size = 1.5, 
              blend = FALSE, 
              interactive = FALSE) + theme(
                panel.background = element_rect(fill = "darkgrey")) 

HelpTDP <-
  DotPlot(object = treatment, 
          features = features,
          cols = c("grey", "#6ab5ba")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(2, 5, 11.5, 12.5), linetype = "dotted", color = "black")+
  geom_vline(xintercept = c(3.75,4.25), linetype = "dotted")

HelpT <- (HelpTLook / HelpTFP) | HelpTDP
HelpT

##### [3.6] Mono & Neutrophils #####
MonoCluster <- c(1,2, 5)
# Create a vector to define colors for each cluster
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[MonoCluster] <- palette.b[MonoCluster]

MonoLook <- 
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

features = c("S100A12", "CXCR1","CCL3","TLR2","BCL2A1", "CCL4L2","NAMPT", "TNFAIP2", "TNFAIP6", "CXCR2", "CXCL8","CCL4", "ITGAX")

MonoFP <- 
  FeaturePlot(object = treatment, 
              features = "TNFAIP2",
              cols = c("lightgrey", "black"),
              label = TRUE,
              pt.size = 1.5, 
              blend = FALSE, 
              interactive = FALSE) + theme(
                panel.background = element_rect(fill = "darkgrey")) 

MonoDP <-
  DotPlot(object = treatment, 
          features = features,
          cols = c("grey", "#900c3e")) +  coord_flip() +  # Flip the x and y a
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), # Adjust x-axis labels angle
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold")) +  # Make y-axis title ("Features") bold) + 
  geom_hline(yintercept = c(0, 2.5, 4.5, 5.5 , 10.5, 11.5 , 13.5, 14.5), linetype = "dotted", color = "black")

Mono <- (MonoLook / MonoFP) | MonoDP
Mono


##### [3.7] Monocytes ####
NCluster <- c(5)
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

##### [4] Altogether #####

one <- 5
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[one] <- palette.b[one]

oneLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )


# monocytes 
two <- c(1, 2)
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[two] <- palette.b[two]

twoLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# B cells 
three <- 7
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[three] <- palette.b[three]

threeLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# CD4 T Helper cells 
four <- 3
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[four] <- palette.b[four]

fourLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# T reg 

# neutrophils
five <- 12
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[five] <- palette.b[five]

fiveLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# CD4 Naive T cells
six <- 4
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[six] <- palette.b[six]

sixLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# CD 8 
seven <- c(6, 8, 9)
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[seven] <- palette.b[seven]

sevenLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

# NK cells
eight <- 10
cluster_colors <- rep("grey", length(unique(treatment.df$seurat_clusters)))
cluster_colors[eight] <- palette.b[eight]

eightLook <- 
  ggplot(treatment.df, aes(x, y, colour = as.factor(seurat_clusters))) +
  geom_point(size = 1, show.legend = FALSE) +
  scale_colour_manual(values = cluster_colors) +
  labs(
    x = "",
    y = "",
    color = ""
  ) + 
  theme_void() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_blank(),  # Remove axis labels
    axis.title = element_blank(),  # Remove axis titles
    plot.margin = margin(1.5, 0.5, 0.5, 0.5, "cm"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right", 
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(1, 0, 0, 0))
  )

##### [5] Combined ####

combined_layout <- (
    
      (oneLook | twoLook|threeLook| four) / 
        (fiveLook | sixLook | sevenLook | eightLook)
    
)

##### [6] Labelled Clusters #####
# Extract UMAP coordinatesb and cluster information
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