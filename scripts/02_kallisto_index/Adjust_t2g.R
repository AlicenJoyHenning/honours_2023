# ADJUSTING T2G FOR SASCRiP

# Load necessary libraries :
library(tidyverse)
library("Biostrings")
?readDNAStringSet()

# Set the working directory
setwd("~/DownloadedMtxTest/AlicenTest/2309114/")

# Direct Adjustment to t2g :  #####
# Read in the output of SASCRiP kallisto_bustools_count (1st output)
transcripts <- read_tsv(
  file = "../2309113/Counts_2309113/Count_analysis/transcripts.txt",
  col_names = FALSE
)
t2gFile <- read_tsv(
  file = "../2309113/New_transcripts_to_genes.txt",
  col_names = FALSE
)

gencodeHeader <- readDNAStringSet(
  file = "../2309113/gencode.v43.transcripts.fa.gz"
)

head(names(gencodeHeader))

Gencodev43Transcripts <- names(gencodeHeader)

Gencodev43df <- data.frame(ENSTLong = Gencodev43Transcripts)

# In the transcripts df - make a new column with transcripts only (the ENST in the t2g)
Gencodev43df <- separate(
  Gencodev43df,
  col = ENSTLong,
  sep = "\\|",
  into = c(paste("Y",1:8)),
  remove = F
)


# Select X1 and Y1
Gencodev43df <- select(
  Gencodev43df,
  c(
    ENSTLong,
    `Y 1`
  )
)


# Rename the transcripts column in t2gFile to Y1
colnames(t2gFile) <- c("ENST", "ENSG", "HGNC")
colnames(Gencodev43df) <- c("ENSTLong", "ENST")

# Now we have overlapping column *ENST* - so we can merge
Gencodet2gMerged <- merge(
  t2gFile,
  Gencodev43df,
  by = "ENST"
)

# Now select and adjust to achieve the correct structure
GencodeT2g <- select(
  Gencodet2gMerged,
  c(
    "ENSTLong",
    "ENSG",
    "HGNC"
  )
)

# Now output this file and get ready to run SASCRiP
write_tsv(
  x = GencodeT2g,
  file = "./GencodeT2g.txt",
  col_names = FALSE
)

