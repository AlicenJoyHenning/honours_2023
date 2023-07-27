# SASCRiP single cell RNA sequencing 
Contents: 
[(i) Installing dependencies](#section-1) 
[(ii) Installing SASCRiP](#section-2) 
[(iii) Preparing fastq files for SASCRiP](#section-3)

__Overview of the pipeline functions :__

- **kallisto_bustools_count** :
INPUT: fastq files (R1 & R2)
OUTPUT: Count_analysis_folder > filtered counts >
filtered_counts.barcodes.txt
filtered_counts.genes.txt
filtered_counts.mtx

Takes as input the short sequencing reads from the experiment and orders them to reflect their location in the reference genome of the organsism. Alignment algorithms compare short reads to a known reference genome (or transcriptome for RNA sequencing). In RNA-seq, psuedoalignment is used: this is when instead of mapping the individual reads directly to the reference transcriptome, the reference transcriptome is manipulated to create an *index*. The reads are then aligned to the *index* using Kallisto, through kb-python. Bustools is then used for gene quantification, in other words after ordering the reads to see which genes they are connected to, you need to quantify how many of each gene (transcripts) there are to be able to analyse gene expression. This is performed at a single cell level using the cell barcodes through bustools.

- **seurat_matrix** :
INPUT:transcripts_to_genes.txt, 
OUTPUT: barcodes.tsv.gz, features.tsv.gz, features_gene_names.tsv, matrix.mtx.gz 
 
Takes the input matrix from kallisto and bustools and the index and makes sure that they are in the correct format for seurat, if not will rearrange them.

- **run_cqc** :
INPUT: matrix in the correct format
OUTPUT: images!! 
Runs Seurart on the matrix to assess quality of cells and remove those that are of low quality, decided by the # of genes expressed and mitochondrial genes expressed. 

- **sctransform_normalize** :
INPUT: seurat object
OUTPUT: images!!
Uses the UMI counts from the healthy cells (after filtered out bad) and generates gene expression values & tells the 2000 most highly variable genes.


## (i) Installing dependencies 
In order for the pipeline to be executed on a device, the following packages must be installed:
+ ***Python*** (greater than version 3.7)
``` python --version ``` (in my case it was Python 3.10.11) if not updated, ``` sudo apt  update ``` / ``` sudo apt install python3 ```
+ ***R*** (greater than version 3.6) ``` R --version ``` (in my case it was R version 4.2.3)
+ ***Seurat R package*** (in my case it wasn't installed yet) R package used for the analysis and visualisation of scRNA-seq data providing a suite of functions including those for quality control and normalisation.
  <pre>
    ```R
    # first check is package is installed:
    library(Seurat)

    # if error message appears, install package: 
    install.packages("Seurat")
    library(Seurat)
    ```
  </pre>
  
+ ***Tidyverse R package*** (in my case this was installed already) R package a collection of R packages for data manipulation
 <pre>
    ```R
    library(tidyverse)
   
   # Otherwise:
   
    install.packages("tidyverse")
    library(tidyverse)
    ```
  </pre>
  
+ ***kb-python*** (**kallisto|bustools**) for pseudo-alignment and quantification capabilities of Kallisto and the single-cell processing utilities of bustools
<pre>
  ```Ubuntu command prompt
  # Check the version available 
  kb --version 

  # If not installed yet, first update **pip** 
# pip is a package manager for python (Pip Installs Packages)
sudo apt update
  sudo apt install python3-pip
  python3 -m pip install --upgrade pip
  
# Install **kb-python** using **pip** 
  # sudo python3 -m pip install kb-python
pip install kb-python # this worked 
  
  # To confirm the installation check the version available 
  kb --version 
```
</pre>
-

## (ii) Installing SASCRiP  
Once the dependcies are installed and updated, it is now time to install the pipeline itself 

NOTE: before using Juypter notebook (lab) on the Ubuntu computer, ensure that the version of R being used by the command prompt (and hence also Juypter) is the one where the R packages above were installed. This can be done by finding the path where the correct version of R is stored and executing the command ```export PATH="/usr/bin:$PATH" ``` then check the version by  opening ``` R ```. After confirmation, close the R with the command ``` q() ``` followed by ``` Save workspace image? [y/n/c]: n ```. 

 The pipeline is run in Juypter-lab. From the command prompt, open using the command ```jupyter-lab``` and uses the package manager ```pip``` : 
 
  ```Jupyter-lab 
  # Install package using pip
  pip install SASCRiP 

  # Install additional R packages required for SASCRiP using the function (not necessary if you've already installed the package, 'install_R_packages'
  import SASCRiP 
  from SASCRiP import sascrip_functions 

  sascrip_functions.install_R_packages()

```

## (iii) Preparing fastq files for SASCRiP 

This step was not necessary for the fastq files used in my project as they were obtained from the 10xv3 chemistry, not the 10xv1 chemistry. This function searches for the RA fastq file that contains both the UMI and transcript sequences that are then separated into their own fastq files to be used as input for the next stage of allignment. 
? Still need to make sure from literature that version 3 was used 

```juypter-lab
edit_10xv1_fastq

import SASCRiP

```

Note that the technology used to generate the fastq files are essental where the format of the fastq file (order of transcript, UMI, and barcode) will be different, where **10xv3** fastqs are ordered as : barcode-UMI, transcript.

Use jupyter, can stay in the same .ipynp throughout, to run the functions of the pipeline: 

# (iv) kallisto_bustools_count

```python

# create variables to set parameters : 
list_of_fastqs = [
    "Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/bamtofastq_S1_L007_R1_001.fastq.gz", 
    "Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/bamtofastq_S1_L007_R2_001.fastq.gz"
]
single_cell_technology = "10Xv3"
output_directory_path = "Alicen/test_files"
species_index = "Alicen/kallisto_index.idx"
species_t2g = "Alicen/transcripts_to_genes.txt"

# alternatively for more than one set of fastqs :

list_of_fastqs = "path/to/FastQ_directory"
input_directory = True
read_separator = ["R1", "R2"]


# Run the function using the variables as inputs : 

sascrip_functions.kallisto_bustools_count(
    list_of_fastqs = list_of_fastqs,
    single_cell_technology = single_cell_technology,
    output_directory_path = output_directory_path, 
    species_index = species_index,
    species_t2g = species_t2g
)                                        

```

# (v) Seurat_matrix 
```python
# create variables for parameters :

matrix_file = "Alicen/test_files/Count_analysis/filtered_counts/filtered_counts.mtx"
gene_index = "Alicen/test_files/Count_analysis/filtered_counts/filtered_counts.genes.txt"
barcode_index = "Alicen/test_files/Count_analysis/filtered_counts/filtered_counts.barcodes.txt"
output_directory = "Alicen/test_files/Count_analysis/filtered_counts"
t2g_file = "Alicen/test_files/Count_analysis/filtered_counts/transcripts_to_genes.txt"

# execute the function:
sascrip_functions.seurat_matrix(
    matrix_file = matrix_file,
    gene_index = gene_index,
    barcode_index = barcode_index,
    output_directory = output_directory, 
    t2g_file = t2g_file, 
    add_hgnc = True
)
```

# (vi) run_cqc
```python

# start with a folder as one of the variables : 

input_file_or_folder = "Alicen/test_files/Count_analysis/filtered_counts/seurat_matrix_output/"
sample_ID = "alpha_test"
output_directory_path = "Alicen/test_files/run_cqc_output/"
gene_column = 2

# execute the function :

sascrip_functions.run_cqc(
    input_file_or_folder = input_file_or_folder,
    sample_ID = sample_ID,
    output_directory_path = output_directory_path,
    gene_column = gene_column
)
```

# (vii) sctransform_normalize
```python
# create variables for parameters :

seurat_object = "Alicen/test_files/run_cqc_output/alpha_test_preQC_seurat.rds"
sample_ID = "alpha_test"
output_directory_path = "Alicen/test_files/sctransform_normalize_ouput/"

# execute the function :
sascrip_functions.sctransform_normalize(
    seurat_object = seurat_object,
    sample_ID = sample_ID,
    output_directory_path = output_directory_path, 
)
```










