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
# In my case:
~/Alicen/tools/bamtofastq_linux/bamtofastq_linux ~/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam ~/Alicen/raw_files/my_fastq_files/human_blood_ifnalpha

```
Alternatively, working on a loop for converting the bam files to fastq files. 

```command prompt
# create array containing the paths to the raw bam files \
bam_files=("/home/intel6700/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam" "/Alic>

# Loop for generating fastq files for each bam file: 

# for bam_file in "${bam_files[@]}"; do 
# output_filename=$(basename "$bam_file")
# extracted_filename="${file_name%%.*}"
# new_filename="home/intel6700/Alicen/raw_files/test_fastq_files/${extracted_filename}"
# /home/intel6700/Alicen/tools/bamtofastq_linux/bamtofastq_linux "${bam_files}" "${new_>
# done 


# loop that creates subfolders with appropriate names:
# for bam_file in "${bam_files[@]}"; do 
# output_filename = "${bam_file%.bam}.fastq"
# mkdir home/intel6700/Alicen/raw_files/test_fastq_files/"${output_filename}"
# home/intel6700/Alicen/tools/bamtofastq_linux/bamtofastq_linux "${bam_files}" "${fastq>
# done 
```



