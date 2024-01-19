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
Once reads have been associated with a tanscriptomic location and a cell of origin, the next step is to determine how many reads are associated with each location/gene in each cell in a process known as quantification. This is helped by the presence of short, random sequences that are unique to each read known as unique molecular identifiers (UMIs). A transcript-to-gene file is used during quantification to associate each transcript (marked by a UMI) with a gene alllowing for number of transcripts associated with each gene, for each cell, to be counted or _quantified_. This file essentially serves as a mapping between the transcriptome and the genes, allowing the software to assign expression values to specific genes based on the aligned or pseudoaligned reads.


1.1 Build _alignment index_ <br>
Psuedoalignment requires the reference transcriptome to be deconstructed and reordered to create what is known as an _alignment index_. Although most psuedoalignment tools are able to autogenerate an alignment index, it is best practice to create one's own using the updated version of the human reference transcriptome. In addition to the index, a transcript to genes (t2g) mapping file must be generated. 

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
