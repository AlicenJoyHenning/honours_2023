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


