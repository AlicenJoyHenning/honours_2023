# Bar graph for receptor expression 
# dots are driving me crazy 
library(Seurat)
library(ggplot2)
library(reshape2)

# [1] Unintentional stacked bar graph ####

receptorinfo <- AverageExpression(TreatmentAnnotated, 
                                  assay = "RNA",
                                  features = c("IFNAR1", "IFNAR2", "IL10RB", "IFNLR1"),
                                  group.by = "celltype",
                                  slot = "data",
                                  verbose = TRUE)
receptorinfo <- as.data.frame(receptorinfo)
receptorinfo <- t(receptorinfo) # get receptors as columns 
rownames(receptorinfo) <- sub("^RNA\\.","", rownames(receptorinfo))
celltype <- rownames(receptorinfo)
celltype <- as.data.frame(celltype)
celltype <- data.frame(lapply(celltype, function(x) gsub("\\.", " ", x))) # removes dots and replaces with spaces for cell types ONLY
receptorinfo <- cbind(celltype, receptorinfo)
receptorinfo$celltype <- factor(receptorinfo$celltype, levels = c(
  "unknown", "monocytes", "neutrophils", "DCs","mDCs","pDCs", "platelets",
  "T helper","naive CD4 T","naive CD8 T","cytotoxic T","NKT","Tregs","Tcm","NK",
  "B"))

# stacked
ggplot(receptorinfo, aes(x = celltype, y = IFNAR1, fill = "IFNAR1")) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  geom_bar(aes(y = IFNAR2, fill = "IFNAR2"), position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex") +
  scale_fill_manual(values = c("IFNAR1" = "#ebebeb", "IFNAR2" = "darkgrey")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
         panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )
# IFNAR1
IFNAR1 <- ggplot(receptorinfo, aes(x = celltype, y = IFNAR1, fill = "IFNAR1", color = "IFNAR1")) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("IFNAR1" = "#efefef")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold",  size = 14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )

# IFNAR2
IFNAR2 <- ggplot(receptorinfo, aes(x = celltype, y = IFNAR2, fill = "IFNAR2", color = "IFNAR2")) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("IFNAR2" = "darkgrey")) +
  scale_color_manual(values = c("IFNAR2" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 14)),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold",  size = 14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )


# IFNLR1
significant_bars <- c(6, 16)  # Indices of bars that are significant

IFNLR1 <- ggplot(receptorinfo, aes(x = celltype, y = IFNLR1, fill = "IFNLR1", color = "IFNLR1")) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("IFNLR1" = "#6ab5ba")) +
  scale_color_manual(values = c("IFNLR1" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 12)),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold",  size = 14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )+ 
  annotate("text", x = significant_bars, y = c(0.24, 0.08), label = "*", vjust = -0.5, size = 8) # significant star 



# IL10RB
IL10RB <- ggplot(receptorinfo, aes(x = celltype, y = IL10RB, fill = "IL10RB", color = "IL10RB")) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  labs(x = "", y = "Average Expression", fill = "Receptor complex", color = "Receptor complex") +
  scale_fill_manual(values = c("IL10RB" = "#b0dfee")) +
  scale_color_manual(values = c("IL10RB" = "black")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 12)),
        legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold",  size = 14),
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )

IFNAR1 / IFNAR2 / IL10RB/ IFNLR1


# [2] Intentional grouped bar graph #####

celltype <- c(rep("unknown", 2), 
              rep("monocytes",2), 
              rep("neutrophils",2), 
              rep("DCs", 2), 
              rep("mDCs", 2), 
              rep("pDCs", 2), 
              rep("platelets",2), 
              rep("T helper",2), 
              rep("naive CD4 T",2), 
              rep("naive CD8 T",2), 
              rep("cytotoxic T",2), 
              rep("NKT",2), 
              rep("Tregs",2), 
              rep("Tcm",2), 
              rep("NK",2), 
              rep("B", 2))
receptor <- rep(c("IFNAR1", "IFNAR2"), 16)
values <- c(receptorinfo$IFNAR1, receptorinfo$IFNAR2)
plotdata <- data.frame(celltype, receptor, values)




# grouped 
ggplot(plotdata, aes(fill = receptor, y = values, x = celltype)) +
  geom_bar(position = "dodge", stat = "identity", width = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()  # Remove minor grid lines
  )
