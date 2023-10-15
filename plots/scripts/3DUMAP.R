# 3D UMAP PLOT 
# CREDIT : https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0/blob/master/3D%20UMAP%20Plotting%20v1.3.R

##### [1] install plotly package : ####
BiocManager::install('plotly')
library(plotly)
library(Seurat)

##### [2] Load & Prepare your Seurat Object : ####
TreatmentAnnotated <- read_rds("alice/honours/results/FinalIndex/TreatmentAnnotated.rds")
TreatmentAnnotated@meta.data$cell_type <- TreatmentAnnotated@active.ident


# Rerun UMAP to get three UMAP axes : 
TreatmentAnnotated <- RunUMAP(TreatmentAnnotated,
                            dims = 1:10,
                            n.components = 3L)

# create meta data label for cell type 
TreatmentAnnotated@meta.data$cell_type <- TreatmentAnnotated@active.ident

# view that labels are UMAP_1 and UMAP_2 and UMAP_3 :
Embeddings(object = TreatmentAnnotated, reduction = "umap")

# Prepare a dataframe for cell plotting 
plot.data <- FetchData(object = TreatmentAnnotated,
                       vars = c("UMAP_1", "UMAP_2", "UMAP_3", "treatment"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))
plot.data$label2 <- paste(TreatmentAnnotated@meta.data$cell_type)
plot.data$label3 <- paste(TreatmentAnnotated@meta.data$treatment)

##### [3] Make plot  ####
# define colour palette for clusters: 
palette.b <- c("#15c284", #0 mono
               "#d72554", #1 CD4_helper
               "#7ac745", #2 neutrophils
               "#9b78c2", #3 T
               # "#9b78c2", #4 mono
               "#e18f75", #5 naive CD8_T
               "#a0d9e9", #6 B 
               "#7e549f", #7 NKT
               "#6ab5ba", #8 mDCs
               "#93aff5", #9 cytotoxic CD8 
               "#c674bc", #10 Tregs
               "#81cfff", #11 NK
               "#74d6c3", #12 platelets
               "#69a923", #13 unknown
               "#00a68e", #14 DCs
               "#dd5839", #15 CD4+
               # "white", #16 CD4 +  
               "#55bbaa") #17 pDC

# define colour palette for treatment type: 
palette.c <- c(
               "#a0d9e9", #1 IFN lambda 
               "#c674bc", #2 IFN alpha
               "grey"  #3 untreated
              ) 


fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~treatment, 
               colors = palette.c,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label3, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text")

# Customize layout to remove gridlines and axes labels
fig <- fig %>%
  layout(
    scene = list(
      xaxis = list(showgrid = FALSE, showticklabels = FALSE, showbackground = FALSE),
      yaxis = list(showgrid = FALSE, showticklabels = FALSE, showbackground = FALSE),
      zaxis = list(showgrid = FALSE, showticklabels = FALSE, showbackground = FALSE)
  )
  
  
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(
          showgrid = FALSE,  # Hide gridlines
          showticklabels = FALSE,  # Hide axis labels
          line = list(color = "white"),  # Change axis color
          gridcolor = "white"  # Change gridline color
        ),
        yaxis = list(
          showgrid = FALSE,
          showticklabels = FALSE,
          line = list(color = "white"),
          gridcolor = "white"
        ),
        zaxis = list(
          showgrid = FALSE,
          showticklabels = FALSE,
          line = list(color = "white"),
          gridcolor = "white"
        )
      )
    )

