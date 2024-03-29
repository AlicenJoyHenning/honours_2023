
![Alt text](https://github.com/AlicenJoyHenning/honours_2023/blob/main/images/00_other/1.png)
## **SINGLE CELL RNA SEQUENCING WORKFLOW**  
## 1. Download and assess the data
### _10X BAM to FASTQ converter_
For the dataset used in this analysis, FASTQ files are not directly available. Rather, the BAM files have been made available that can be downloaded and converted into FASTQ files using the flowing steps [(download data)](https://github.com/AlicenJoyHenning/honours_2023/edit/main/docs/00_downloading_files.md) 
with the use of `bamtofastq_linux`. 

Before running the downstream analysis, some quality control checks need to be done to ensure the raw data has no underlying problems, or inform problems that you may have.
`FastQC` provides a QC report, summarised in `MultiQC`, that can be run in a non-interactive mode where it would be suitable for integrating into a larger analysis pipeline [(guide](https://github.com/AlicenJoyHenning/honours_2023/blob/main/docs/01_sequencing_quality_assessment.md);[ & code here)](scripts/01_seqQC/quality_control.sh).
> _STUDY OUTCOME 1 | Initial quality assessment_<br>
>
> Across the investigated datasets, the quality of sequencing data was deemed acceptable by FastQC metrics.
> Although a subset of sequences failed to reach certain quality thresholds set by FastQC, this could be attributed to the nature of the 10X Genomics Chromium scRNA-seq output.
> Specifically, an output known as the index file contains intentionally identical sequences to allow the Illumina sequencing technology to differentiate between adjacent read pairs.
> Consequently, these sequences are highly over-represented in the read files and exhibit a non-normal distribution of bases and GC content, as was detected by FastQC.
> However, a portion of the failed sequences, less than 15 %, arose from the true sequencing reads.
> This was observed in the IFN-treated samples and, while not a concern worthy of disregarding the datasets, was kept in mind throughout downstream analysis.
><br>
<br>

## 2. Pre-psuedoalignment processes
### _Building the reference transcriptome input for *Kallisto* to be used in the SASCRiP pipeline_
#### (i) Install and download the necessary dependencies
To run the pseudo alignment tool (**Kallisto**), an index of the reference transcriptome is needed. Although the SASCRiP function **kallisto_bustools** is able to do this automatically by changing some parameters, I needed to know how to do this manually.<br> 
For this process, python needs to be installed along with the package manager, pip. Download the latest version of python (<a href="https://www.python.org/downloads/" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">PYTHON</a>) and download the file here: (<a href="https://bootstrap.pypa.io/get-pip.py" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">PIP</a>). Ensure Python is installed from the cmd using ```python --version``` before installing pip ```python get-pip.py``` which can be verified afterwards as well using ```pip --version```. Next, JupyterLab must be installed which can now be done from the cmd ```pip install jupyterlab``` and opened  ```jupyter lab```. Note that this takes you to a Google Chrome page with Jupyter ready to be used, this is where the second portion of this tutorial needs to be completed.

This process also requires Kallisto to be installed (<a href="https://github.com/pachterlab/kallisto/releases" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">Kallisto</a>). From the options, choose and install the one compatible with your device and to use it, travel to the directory where it is saved ```cd kallisto_windows-v0-50.0\kallisto``` : 
+ kallisto_linux-v0.50.0.tar.gz
+ kallisto_mac-v0.50.0.tar.gz
+ kallisto_mac_m1-v0.50.0.tar.gz
+ kallisto_windows-v0.50.0.zip<br>

From this point, the full transcriptome from Ensembl (files ending in cdna.all.fa.gz) must be downloaded. To build the human transcriptome index, first download the transcriptome, which is available under cDNA on the [Ensembl website](http://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/) 
and execute [the following](https://github.com/AlicenJoyHenning/honours_2023/blob/main/scripts/02_kallisto_index/download_transcriptome.sh) in the command prompt. 

> _STUDY OUTCOME 2 | Successful download of the full transcriptome_<br>
> <div style="text-align:center">
>    <img src="https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/7c5f8b9b-e275-4dd4-b79f-1ddc0e55b37f" alt="image">
> </div>
> NOTE: you could, of course, do this by clicking the download button when you travel to the website.<br>
<br>

#### (ii) Build the index file 
After downloading the full transcriptome file, you now need to build the index itself. This was done in the command prompt using Kallisto [as follows](https://github.com/AlicenJoyHenning/honours_2023/blob/main/scripts/02_kallisto_index/kallisto_index.sh).

> _STUDY OUTCOME 3 | Successful building of the pseudoalignment index_<br>
> ![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/ecef3b27-fc09-4847-801e-42bc643877d6)<br>


#### (iii) Build the transcripts to gene file 
Once the index is created, the transcripts to genes text file must also be compiled. This can be done using a function from `kb_python` called `create_t2g_from_gtf` . This requires gtf (gene transfer format) files as input that must be downloaded. To download gtf files go to ensembl website > human > latest genome assembly > GRCh38 (or latest version) > access the gtf file (Homo_sapiens.GRCh38.110.gtf.gz). Once downloaded, store the gtf file in a specific directory and complete [the following](scripts/02_kallisto_index/manual_kb_python.ipynb). 

#### (iv) Adjust the transcripts to gene file for Kallisto 


### 3. SASCRiP: psuedoalignment and quantification pipeline 

### 4. Cell and data quality control 

### 5. Integration of study samples

### 6. Dimensionality reduction and clustering of integrated dataset

### 7. Differential gene expression analyses 
