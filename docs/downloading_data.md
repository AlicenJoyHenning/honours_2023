# Downloading the data

## 10X BAM to FASTQ converter

### Download BAM files 

Download the bam files from the ena database, create a text file with the links to the bam files :  
Untreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769384/human_blood_untreated.bam.bam  
IFNalphatreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769382/human_blood_ifnalpha.bam.bam  
IFNlambdatreated: ftp://ftp.sra.ebi.ac.uk/vol1/run/SRR107/SRR10769383/human_blood_ifnlambda.bam.bam

```command prompt
# Travel to the directory you want the bam files to be stored (alternatively you could enter the path manually) : 
cd Alicen/raw_files/bam_files

# Create a text file containing the urls for the files to be downloaded with each url on a new line :
nano bam_urls.txt

# Use wget to download the bam files from the list using the -i option : 
wget -i bam_urls.txt

```
### Convert BAM files to FASTQ files  

Install the bamtofastq_linux (https://github.com/10XGenomics/bamtofastq) and store it in a directory of your choice. After installing bamtofastq_linux (should appear with a gear symbol), travel to the directory in the command prompt and execute the following : 

```command prompt
# Travel to the correct directory : 
cd /path/to/bamtofastq_linux

# To convert the downloaded bam_files into fastq files one by one :
/path/to/bamtofastq_linux /path/to/bam_files /path/to/fastq_files

# In my case :
~/Alicen/tools/bamtofastq_linux/bamtofastq_linux ~/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam ~/Alicen/raw_files/my_fastq_files/human_blood_ifnalpha

```
Alternatively, this loop converts the bam files to fastq files without having to enter in the paths to each bam file individually : 

```Linux command prompt
# Create array containing the paths to the raw bam files :
bam_files=("/home/intel6700/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam" "/home/intel6700/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam" "/home/intel6700/Alicen/raw_files/bam_files/human_blood_untreated.bam")


# Loop for generating fastq files for each bam file: 
for bam_file in "${bam_files[@]}"; do 
output_filename = basename "$bam_file"%.bam
/home/intel6700/Alicen/tools/bamtofastq_linux/bamtofastq_linux "${bam_files}" ~/home/intel6700/Alicen/raw_files/fastq_files/"${output_filename}"
done 
