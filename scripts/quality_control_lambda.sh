
# Loop that allows you to enter a single command as input for fastqc analysis instead of manually inputting individual fastq files 

fastq_directory="/home/intel6700/Alicen/raw_files/fastq/human_blood_ifnlambda/human_blood_ifnlambda"
output_directory="/home/intel6700/Alicen/quality_control/ifnlambda"

for fastq_file in ${fastq_directory}/*.gz; do
fastqc $fastq_file -o $output_directory
done 




