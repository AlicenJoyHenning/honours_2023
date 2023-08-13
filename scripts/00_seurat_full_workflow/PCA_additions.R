# creating a new section in the seurat object to manually store the PCA data : 
# Create the data frame
alpha_positive_genes_PC1 <- c(
  "RPL23A", "RPL3", "RPL13A", "RPL26", "RPS3A",  "RPL18A", "RPS23", "RPL7A", "RPS18", "RPL35",
  "ENSG00000237550", "RPSA", "RPL10A", "RPL5", "RPL10",  "RPS4X", "EEF1B2", "RPL6", "RPS7", "RPS8",
  "RPS5", "RPL18", "RPS20", "RPL37A", "RPS15A", "RPL12", "RPL24", "RPL35A", "RPS6", "RPS3"
)

alpha_negative_genes_PC1 <- c(
  "IL1RN", "S100A8", "FTH1", "FFAR2", "MNDA", "ACSL1", "APOBEC3A", "TYROBP", "NCF1C", "NINJ1","GLUL", "SNX10", "IFITM3", "FCER1G", "SPI1", 
  "NCF1", "TYMP", "C15orf48", "IER3", "TLR2", "CCL4", "CLEC7A", "CCL4L2", "LST1", "S100A11", "NCF2", "VNN2", "SAMSN1", "ISG15", "BASP1"
)

alpha_positive_genes_PC2 <- c("LEF1", "LDHB", "IL7R", "BCL11B", "IKZF1", "TCF7", "GPR155", "ADTRP", "ITGA6", "CAMK4",
                              "GIMAP4", "CMTM8", "LEPROTL1", "C12orf57", "MAL", "PRMT2", "PDCD4", "CD3G", "CCR7", "TNFRSF25", 
                              "MYC", "SCML1", "LTB", "STING1", "GIMAP5", "FHIT", "FAM184A", "GIMAP6", "PASK", "ABLIM1") 


alpha_negative_genes_PC2 <- c("HLA-DPA1", "HLA-DRA", "CD74", "MYOF", "KYNU", "HLA-DQB1", "LILRB1", "CD86", "RGL1", "HLA-DQA1",
                              "LGALS1", "ADAM28", "EPB41L3", "CPVL", "CYP1B1", "OLR1", "ANXA2", "LACC1", "EMP1", "LY86", 
                              "BANK1", "MS4A1", "TCF4", "MS4A7", "RIN2", "CCL2", "CCL7", "SWAP70", "ENG", "FGD2")  

alpha_positive_genes_PC3 <- c("MS4A1", "IGHM", "BANK1", "BCL11A", "IGHD", "IGKC", "HLA-DQA1", "HLA-DPA1", "ADAM28", "HLA-DRA", 
                              "TCF4", "IGKJ5", "CCSER1", "CD19", "PLEKHG7", "RAB30", "BCL2", "MEF2C", "NAPSB", "EAF2", 
                              "HLA-DQB1", "CCR7", "IGLC2", "CD79B", "HLA-DPB1", "STAP1", "CD24", "ACSM3", "CD74", "C7orf50")

alpha_negative_genes_PC3 <- c("NKG7", "PRF1", "GNLY", "ANXA1", "GZMB", "CLIC3", "SAMD3", "MYBL1", "METRNL", "SYTL2", 
                              "CD8A", "FYN", "SH2D2A", "C1orf21", "KLRD1", "IL18RAP", "GZMM", "CD96", "PYHIN1", "PLAAT4", 
                              "SYTL3", "AUTS2", "EOMES", "KLRF1", "IQGAP2", "PIK3R1", "LAG3", "MAF", "TNFRSF18", "TRGC2")
alpha_positive_genes_PC4 <- c("NKG7", "PRF1", "GNLY", "GZMB", "CLIC3", "SAMD3", "KLRD1", "MS4A1", "KLRF1", "BANK1", 
                              "IGHM", "BCL11A", "IL18RAP", "MYBL1", "CD38", "C1orf21", "IGHD", "SH2D2A", "EOMES", "ADAM28", 
                              "CD8A", "IGKC", "STAP1", "TNFRSF18", "SYTL2", "CCSER1", "PYHIN1", "FGFBP2", "GZMM", "AUTS2") 

alpha_negative_genes_PC4 <- c("MYOF", "EPB41L3", "RGL1", "CYP1B1", "CPVL", "EMP1", "OLR1", "CCL2", "CCL7", "LACC1", 
                              "CD86", "CXCL11", "ENG", "MS4A7", "PLA2G7", "CTSB", "P2RY6", "LYZ", "MSR1", "LILRB4", 
                              "TGFBI", "LILRB1", "LEF1", "SRC", "LGALS1", "CST3", "KYNU", "RIN2", "THBS1", "RGS10")

alpha_positive_genes_PC5 <- c("ENSG00000284874", "NRGN", "ACRBP", "TSC22D1", "DAB2", "PTCRA", "TUBA4A", "CMTM5", "MMD", "TMEM40", 
                              "C2orf88", "ENKUR", "CLU", "PGRMC1", "CTTN", "H2BC11", "MPIG6B", "MTURN", "SPARC", "PDE5A", 
                              "GTF3C2", "MAP3K7CL", "H2AJ", "TNFSF4", "CLDN5", "BEX3", "H2AC6", "CTSA", "SSX2IP", "TSPAN33") 

alpha_negative_genes_PC5 <- c("ISG15", "MYOF", "RGL1", "LGALS1", "S100A11", "APOBEC3A", "LILRB1", "TYROBP", "OLR1", "MS4A7", 
                              "LACC1", "RIN2", "LILRB4", "CPVL", "KYNU", "EMP1", "IL1RN", "ANXA5", "CTSB", "S100A4", 
                              "IFITM3", "ENG", "CYP1B1", "MT2A", "CCL4", "PLA2G7", "S100A6", "GCH1", "IL4I1", "CD86") 

alpha_pca_df <- data.frame(
  PCA = rep(1:5, each = 60),
  Feature = c(alpha_positive_genes_PC1, alpha_negative_genes_PC1, alpha_positive_genes_PC2, alpha_negative_genes_PC2,alpha_positive_genes_PC3, alpha_negative_genes_PC3, alpha_positive_genes_PC4, alpha_negative_genes_PC4, alpha_positive_genes_PC5, alpha_positive_genes_PC5),
  "Pos/Neg" = rep.int(rep(c("Positive", "Negative"), each = 30),5)
)


alpha@assays$PCA <- alpha_pca_df