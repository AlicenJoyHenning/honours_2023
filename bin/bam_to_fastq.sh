# Using samtools, downloaded BAM files are converted into fastq files 

# First check the quality of each bam file:
samtools quickcheck human_blood_ifnalpha.bam 
samtools quickcheck human_blood_ifnlambda.bam
samtools quickcheck human_blood_untreated.bam 

# Loop that takes as input the path to a bam file and gives as output two paired end fastq files : 
# The assumes you are in the directory where the bam files are located 
bam_files = ("human_blood_ifnalpha.bam" "human_blood_ifnlambda.bam" "human_blood_untreated.bam")
output_directory = "Alicen/raw_files/fastq_files"

///

# Loop for generating 1 fastq file for each bam file: 

# for bam_file in "${bam_files[@]}"; do 
# output_filename = "${bam_file%.bam}.fastq"
# mkdir ~/Alicen/raw_files/fastq_files/"${output_filename}"
# samtools fastq -@ 4 "${bam_file}" > "${fastq_output}/${output_filename}"
# done 

///

# Loop for generating 2 paired end fastq files for each bam file (needed as input for SASCRiP pipeline): 
for bam_file in "${bam_file[@]}"; do 
output_filename = "${bam_file%.bam}.fastq"
mkdir ~/Alicen/raw_files/fastq_files/"${output_filename}"
samtools fastq -@ 4 -1 
