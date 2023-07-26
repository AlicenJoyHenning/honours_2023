# Downloading the data

## 10X BAM to FASTQ converter

Download the bam files from the ena database, create a text file with the links to the bam files: 
Untreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769384/human_blood_untreated.bam.bam
IFNalphatreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769382/human_blood_ifnalpha.bam.bam
IFNlambdatreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769383/human_blood_ifnlambda.bam.bam

```command prompt
nano bam_urls.txt 
wget <url>
```

Install the bamtofastq_linux and set it in a certain directory. After installing the bamtofastq_linux, travel to the directory. 

```command prompt
cd /path/to/bamtofastq_linux

# To convert the downloaded bam_files into fastq files one by one:
/path/to/bamtofastq_linux /path/to/bam_files /path/to/fastq_files

# Alternatively, you can use a loop to execute them simultaneously:
bam_files = [
for 
```

