# Kallisto Index: transcript to genes file with updated information 

# Load dependencies 
BiocManager::install('BUSpaRse')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
library(BUSpaRse)
library(BSgenome.Hsapiens.UCSC.hg38)

# Use the tr2g_ensembl function to create the t2g file 
?tr2g_ensembl()

tr2g_hs <- tr2g_ensembl("Homo sapiens", 
                        use_gene_version = FALSE, 
                        gene_biotype_use = "cellranger", 
                        write_tr2g = TRUE,
                        out_path = "honours/kallisto_index/01.09.version/")
# 