# B Cell differentiation markers 

library(Seurat)
library(ggplot2)
library(reshape2)

# [0] Calculate expression 
# isolate B cells
Bcells <- subset(TreatmentAnnotated, celltype == "B")

markerinfo <- AverageExpression(Bcells, 
                                  assay = "RNA",
                                  features = c("OAS2", "IFIT3", "STAT1", "CD38", "JCHAIN", "XBP1"),
                                  group.by = "sample",
                                  slot = "data",
                                  verbose = TRUE)


markerinfo <- as.data.frame(markerinfo)
markerinfo <- t(markerinfo) # get receptors as columns 
rownames(markerinfo) <- sub("^RNA\\.","", rownames(markerinfo))
sample <- rownames(markerinfo)
sample <- as.data.frame(sample)
sample <- data.frame(lapply(sample, function(x) gsub("\\.", " ", x))) # removes dots and replaces with spaces for cell types ONLY

markerinfo <- cbind(sample, markerinfo)

# [1] 

OAS2 <- ggplot(markerinfo, aes(x = sample, y = OAS2, fill = "OAS2", color = "OAS2")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("OAS2" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  ggtitle("OAS2") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

IFIT3 <- ggplot(markerinfo, aes(x = sample, y = IFIT3, fill = "IFIT3", color = "IFIT3")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("IFIT3" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  ggtitle("IFIT3") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

STAT1 <- ggplot(markerinfo, aes(x = sample, y = STAT1, fill = "STAT1", color = "STAT1")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("STAT1" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  ggtitle("STAT1") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))



OAS2 | IFIT3 | STAT1

CD38 <- ggplot(markerinfo, aes(x = sample, y = CD38, fill = "CD38", color = "CD38")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("CD38" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  ggtitle("CD38") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

JCHAIN <- ggplot(markerinfo, aes(x = sample, y = JCHAIN, fill = "JCHAIN", color = "JCHAIN")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("JCHAIN" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  ggtitle("JCHAIN") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

XBP1 <- ggplot(markerinfo, aes(x = sample, y = XBP1, fill = "XBP1", color = "XBP1")) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("XBP1" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(face = "bold", size = 14), # angle = 45, hjust = 1, 
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.position = "none",
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+
  ggtitle("XBP1") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

CD38 | JCHAIN | XBP1
