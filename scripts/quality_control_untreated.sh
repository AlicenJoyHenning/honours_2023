
# Loop that allows you to enter a single command as input for fastqc analysis instead of manually inputting individual fastq files 

fastq_directory="/home/intel6700/Alicen/raw_files/fastq/human_blood_untreated/human_blood_untreated"
output_directory="/home/intel6700/Alicen/quality_control/untreated"

for fastq_file in ${fastq_directory}/*.gz; do
fastqc $fastq_file -o $output_directory
done 




