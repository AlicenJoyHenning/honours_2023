# FastQC: Sequencing Quality Control

Before running the downstream analysis, some quality control checks need to be done to ensure the raw data has no underlying problems, or inform problems that you may have.
**FastQC** provides a QC report that can be run in a non-interactive mode where it would be suitable for integrating into a larger analysis pipeline. 

```Unix Command Prompt
# check is the FastQC is installed on the system : 
fastqc

# If not, then install it using the following (for an Ubuntu system) : 
sudo apt uupdate 
sudo apt install fastqc

# Once FastQC is installed on the Ubuntu system, it can be used to analyse the quality of FASTQC files using the comand :
fastqc /path/to/file.fastq
```

# MultiQC : collecting results 

```linux
pip install multiqc

# check :
multiqc --version
# output : 1.14

```

Create a script to run multiQC using input from fastQC outputs :

``` commandline
# Basic Premise: 
multiqc input/directoy -o output/directory 
multiqc quality_control/ifnalpha/fastQC -o quality_control/ifnalpha/multiQC
```
Analyze the output : 
I1 : these reads contain index sequences (cell barcodes) that are used to assign the individual cells sequenced in the experiment 
R1: sequence data from cDNA fragments originating from the 5' end 
R2 : sequence data from cDNA fragments originating from the 3' end 
