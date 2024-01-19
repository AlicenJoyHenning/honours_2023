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
scRNAseq FASTQ files hold the raw sequencing data (reads) that on their own hold no biological insight. The reads must be aligned to the reference genome to _piece together_ the picture of the transcritpome they represent. Alignment is approximated by two techniques that differ in both accuracy and computational resources. The more computationally intensive genome alignment technique aligns reads across the full length of the reference genome. Alternatively, psuedoalignment deconstructs the reference genome and aligns reads to smaller segments in an optimized order making the process more _compuatationally digestable_ albeit at a decrease in the accuracy of alignment. In scRNAseq, this decrease is not detrimental to downstream analysis and, considering the overwhelming computational resources of scRNAseq data when compared to bulk RNAseq data, makes scRNAseq analyses more accessible to the scientific community at large. 

1.1 Build alignment index<br>

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
