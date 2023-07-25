# To convert files from an sra file format into paired-end fastq files

# Install the sratoolkit : 
# GitHub - s-andrews/sradownloader: A script to make downloading of SRA/GEO data easier 
# Travel to the directory in which it is located : 
# cd Alicen\tools\sratoolkit\sratoolkit.v1.1.2\bin
# Configure the toolkit: 
./vdb-config -i

# Command to download the sra files using their accession numbers :  
# prefetch --output-directory /path/to/output/directory accession_number
prefetch --output-directory ~/Alicen/raw_files/sra_files SRR10769384

# Once downloaded and stored in a directory, the sra files were converted to fastq files: 
fasterq-dump SRR10769384 --split-3 -o ~/Alicen/raw_files/fastq_files 

