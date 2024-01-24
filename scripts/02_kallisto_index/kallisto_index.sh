# Run kallisto index. kallisto will work on .fa and .fz.gz files so there is no need to unzip the downloaded file:

kallisto index -i 	Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx	Homo_sapiens.GRCh38.cdna.all.fa.gz
# -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx: specifies the output index file name
# Homo_sapiens.GRCh38.cdna.all.fa.gz: input transcriptome reference file where -i flag indicates that you want to build an index from this file
# OUTPUT :
[build] building MPHF
[build] creating equivalence classes ...
[build] target de Bruijn graph has k-mer length 31 and minimizer length 23
[build] target de Bruijn graph has 936536 contigs and contains 108619921 k-mers 
