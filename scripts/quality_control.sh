
# Loop that allows you to enter a single command as input for fastqc analysis instead of manually inputting individual fastq files 

fastq_directory = "/home/intel6700/Alicen/raw_files/fastq/human_blood_ifnalpha"
output_directory = "/home/intel6700/Alicen/quality_control/ifnalpha"

for fastq_file in "$fastq_directory"; do
fastqc "$fastq_file" -o "$output_directory"
done 


# This is able to find all the fastq files and store them in a particular directory : 
# find . -name "*fastq.gz"> ~/Alicen/quality_control/ifnalpha/output.txt

# To run a check on one file : 
# fastqc name.fastq.gz -o ~/output/directory

