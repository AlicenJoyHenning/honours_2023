## Overview of scRNAseq analysis 
[**Contents:**](#section-1) <br>
_1. Preprocessing<br>
&nbsp;&nbsp;&nbsp;1.1 Build alignment index<br> 
&nbsp;&nbsp;&nbsp;1.2 Psueoalign<br>
&nbsp;&nbsp;&nbsp;1.3 Quantification
<br>&nbsp;&nbsp;&nbsp;1.4 Remove low quality cells
<br>&nbsp;&nbsp;&nbsp;1.5 Normalise
<br>&nbsp;&nbsp;&nbsp;1.6 Scale data
<br>&nbsp;&nbsp;&nbsp;1.7 Integrate<br>
2. Identify and visualise cell populations<br>
&nbsp;&nbsp;&nbsp;2.1 Linear dimensionality reduction (PCA)
<br>&nbsp;&nbsp;&nbsp;2.2 Clustering (KNN)
<br>&nbsp;&nbsp;&nbsp;2.3 Uniform Manifold Approximation and Projection (UMAP)<br>
3. Annotate cell populations<br>
4. Investigate gene expression of cell populations<br>_

[**1. Preprocessing**](#section-1)<br>
READ > TRANSCRIPTOME LOCATION<br>
scRNAseq FASTQ files hold the raw sequencing data (reads) that on their own hold no biological insight. The reads must be aligned to the reference genome to _piece together_ the picture of the transcritpome they represent. Alignment is approximated by two techniques that differ in both accuracy and computational resources. The more computationally intensive genome alignment technique aligns reads across the full length of the reference genome. Alternatively, psuedoalignment deconstructs the reference genome and aligns reads to smaller segments in an optimized order making the process more _compuatationally digestable_ albeit at a decrease in the alignment rate. In scRNAseq, this decrease is not detrimental to downstream analysis and, considering the overwhelming computational resources of scRNAseq data when compared to bulk RNAseq data, makes scRNAseq analyses more accessible to the scientific community at large. For these reasons, this analysis utilises psuedoalignment rather than genome alignment.

READ + TRANSCRIPTOMIC LOCATION > CELL OF ORIGIN<br>
In single cell RNA sequencing, reads must be reassigned to the single cell from which they came. According to the technology used, short sequences are incorporated into sequencing reads that originated from the same cell (10XGenomics use cell barcodes). Using these sequences, reads from the same cell, having the same identifier, are grouped together.

HOW MANY READS IN EACH TRANSCRIPTOMIC LOCATION FOR EACH CELL<br>
_COUNT MATRIX_<br>
Once reads have been associated with a tanscriptomic location and a cell of origin, the next step is to determine how many reads are associated with each location/gene in each cell in a process known as quantification. This is helped by the presence of short, random sequences that are unique to each read known as unique molecular identifiers (UMIs). A transcript-to-gene file is used during quantification to associate each transcript (marked by a UMI) with a gene alllowing for number of transcripts associated with each gene, for each cell, to be counted or _quantified_. This file essentially serves as a mapping between the transcriptome and the genes, allowing the software to assign expression values to specific genes based on the aligned or pseudoaligned reads. Ultimately, the data is outputted in the form of a _cell x gene_ matrix. 

Knowing this preliminary gene expression data for each cell allows decisions to be made surrounding quality control. In scRNAseq, cells can contain low quality data for a number of reasons, some of which are detectable and targeted for quality control such as mitochondrial percentage. Generally speaking healthy cells can be differentiated by their intact membranes. However, the integrity of the cell membrane of dead and dying cells is often comprimised meaning the RNA present in the cytoplasm may escape. Comparitively, the RNA present in the mitochondria remains fairly consistent. The ratio of cytoplamsic RNA to mitochondrial RNA will therefore change in accordance to cell health where unhealthy/dead cells likely have higher percentage of mitochondrial RNA (above 10 %). In scRNAseq, controls also need to be taken against events that may occur during sequencing such as the sequencing of two cells as one (doublets) or the sequencing of nothing. Although an ongoing problem, these cases can be approcimately identified using a statistical method of outlier identification. These may include standard deviation (SD) or median absolute deviation (MAD), which is more resilient to outliers in a dataset (SD squares distances from the mean so large deviations are weighted more heavily disguising true otuliers). 




1.1 Build _alignment index_ <br>
Psuedoalignment requires the reference transcriptome to be deconstructed and reordered to create what is known as an _alignment index_. Although most psuedoalignment tools are able to autogenerate an alignment index, it is best practice to create one's own using the updated version of the human reference transcriptome.  

1.2 Psuedoalign<br>

1.3 Quantification<br>

1.4 Remove low quality cells<br>

1.5 Normalise<br>

1.6 Scale data<br>

1.7 Integrate<br>

[**2. Identify and visualise cell populations**](#section-1)<br>
2.1 Linear dimensionality reduction (PCA)<br>

2.2 Clustering (KNN)<br>

2.3 Uniform Manifold Approximation and Projection (UMAP)<br>

[**3. Annotate cell populations**](#section-1)<br>

[**4. Investigate gene expression of cell populations**](#section-1)<br>
