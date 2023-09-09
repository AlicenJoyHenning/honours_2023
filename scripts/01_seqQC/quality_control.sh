
# SEQUENCING QUALITY CONTROL 

# To check if fastqc and multiqc are functional use one file : 
# fastqc name.fastq.gz -o ~/fastqc/output/directory
# multiqc . -o ~/multiqc/output/directory

# Loop that allows you to enter a single command as input for fastqc analysis instead of manually inputting individual fastq files 

fastq_directory="/home/intel6700/Alicen/raw_files/fastq/human_blood_ifnalpha"
fastqc_output_directory="/home/intel6700/Alicen/quality_control/ifnalpha/fastqc"
mutliqc_output_directory="/home/intel6700/Alicen/quality_control/ifnalpha/multiqc"

for fastq_file in ${fastq_directory}/*.gz; do
fastqc $fastq_file -o $fastqc_output_directory
multiqc $fastqc_output_directory -o $multiqc_output_directory
done 

echo "FastQC and MultiQC analysis completed."




