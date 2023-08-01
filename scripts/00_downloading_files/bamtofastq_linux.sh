# create array containing the paths to the raw bam files \
bam_files=("/home/intel6700/Alicen/raw_files/bam_files/human_blood_ifnalpha.bam" "/Alicen/raw_files/bam_files/human_blood_ifnlambda.bam" "/Alicen/raw_files/bam_files/human_blood_untreated.bam") 

# Loop for generating fastq files for each bam file: 

# for bam_file in "${bam_files[@]}"; do 
# output_filename=$(basename "$bam_file")
# extracted_filename="${file_name%%.*}"
# new_filename="home/intel6700/Alicen/raw_files/test_fastq_files/${extracted_filename}"
# /home/intel6700/Alicen/tools/bamtofastq_linux/bamtofastq_linux "${bam_files}" "${new_filename}"
# done 


# loop that creates subfolders with appropriate names:
# for bam_file in "${bam_files[@]}"; do 
# output_filename = "${bam_file%.bam}.fastq"
# mkdir home/intel6700/Alicen/raw_files/test_fastq_files/"${output_filename}"
# home/intel6700/Alicen/tools/bamtofastq_linux/bamtofastq_linux "${bam_files}" "${fastq_output}" 
# done 
