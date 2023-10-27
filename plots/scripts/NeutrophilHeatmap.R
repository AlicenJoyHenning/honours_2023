# Neutrophil  markers 

library(Seurat)
library(ggplot2)
library(reshape2)

# [0] Calculate expression 
# isolate B cells
Neutrophils <- subset(TreatmentAnnotated, celltype == "neutrophils")

markerinfo <- AverageExpression(Neutrophils, 
                                assay = "RNA",
                                features = c(
                                             "OAS1", "OAS2",
                                            "IFIT1", "IFIT2",  
                                            "STAT1", "STAT2",
                                            "MX1",
                                            "THBD","CAMP", "CD93", "VDR",
                                            "IFNAR1", "IFNLR1"
                                            ),
                                group.by = "treatment",
                                slot = "data",
                                verbose = TRUE)


markerinfo <- as.data.frame(markerinfo)
markerinfo <- t(markerinfo) # get receptors as columns 
rownames(markerinfo) <- sub("^RNA\\.","", rownames(markerinfo))
sample <- rownames(markerinfo)
sample <- as.data.frame(sample)
sample <- data.frame(lapply(sample, function(x) gsub("\\.", " ", x))) # removes dots and replaces with spaces for cell types ONLY

markerinfo <- cbind(sample, markerinfo)

# Heatmap 

hot <- data.frame(
  sample = c("α", "λ", "u"),
  IFNAR1 = c(35, 40, 45),
  MX1 = c(38, 10, 0),
  OAS1 =  c(30, 0, 0),
  IFIT1 =  c(42, 0, 0),
  STAT1 = c(25,  0,  0),
  STAT3 = c(30,  0,  0),
  
  IFNLR1 = c(10, 20, 15),
  VDR = c(0,  35,  20),
  CD93 = c(0, 38, 34),
  THBD = c(0,  45,  0)
)

# Reorder the GO column 
hot$sample <- factor(hot$sample, levels = c(c("u","λ","α")))

hotter <- melt(hot, id.vars = "sample", variable.name = "Gene", value.name = "Count")


ggplot(hotter, aes(x = Gene, y = sample, fill = Count)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "#a298bd") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5),  # Rotate x-axis labels and set font size and style
        axis.text.y = element_text(size = 14, angle = 90, vjust = 0.5),
        axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
        legend.key.size = unit(2, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = "bold")) +
  labs(x = "", y = "", fill = "Count") +
  scale_x_discrete(labels = c("IFNAR1", "MX1", "OAS1", "IFIT1", "STAT1", "STAT3", "IFNLR1", "VDR", "CD93", "THBD")) +
  theme(legend.position = "left")
