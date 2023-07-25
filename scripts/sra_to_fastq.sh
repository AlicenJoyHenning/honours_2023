# To convert files from an sra file format into paired-end fastq files

# Install the sratoolkit : 
# GitHub - s-andrews/sradownloader: A script to make downloading of SRA/GEO data easier 
# Travel to the directory in which it is located : 
# cd Alicen\tools\sratoolkit\sratoolkit.v1.1.2\bin
# Configure the toolkit: 
./vdb-config -i


////

# Command to download the sra files individually using their accession numbers :  
# prefetch --output-directory /path/to/output/directory accession_number
# prefetch --output-directory ~/Alicen/raw_files/sra_files SRR10769384

# Once downloaded and stored in an output directory, the sra files were converted to two fastq paired end reads files
# input_reads_1.fastq (forward) and input_reads_2.fastq (reverse)
# fasterq-dump SRR10769384 --split-3 -o ~/Alicen/raw_files/fastq_files 

////

# Loop to download the sra files altogether and convert them to paired end fastq files using their accession numbers :
# Note that both output directories must exist already 

sra_input = ("SRR10769382" "SRR10769383" "SRR10769384")  # defining the SRA accession numbers as an array
sra_output_directory = "Alicen/raw_files/sra_files" # defining the output directories 
fastq_output_directory = "Alicen/raw_files/fastq_files"

for sra_accession in "${sra_input[@]}; do 
prefetch --output-directory "${sra_output_directory}" "${sra_accesion}" 
fasterq-dump  "${sra_output_directory}/${sra_accession}.sra" --split-3 -o "${fastq_output_directory}"
done


////

# Suggestion loop: 
# Loop over each accession number
# 
# for accession in $accession_numbers; do
#    # Download SRA file
#    prefetch $accession
#    
#    # Convert SRA file to FASTQ files
#    fastq-dump --split-files $accession
#    
#    # Optional: Remove the SRA file to save disk space
#    rm ${accession}.sra
# done
