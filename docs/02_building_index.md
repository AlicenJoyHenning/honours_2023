
# Building the reference transcriptome input for kallisto 

To run the pseudo alignment tool (kallisto), an index of the reference transcriptome is needed. Although the SASCRiP function **kallisto_bustools** is able to do this automatically by changing some parameters, I needed to know how to do this manually: 

The full transcriptome from Ensembl (files ending in cdna.all.fa.gz) must be downloaded. To build the human transcriptome index, first download the transcriptome, which is available under cDNA on the Ensembl website, at http://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/, and execute the following in the command prompt : 

  
```command promt 
# Download the full transcriptome from ensemble : 
curl -O ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Run kallisto index. kallisto will work on .fa and .fz.gz files so there is no need to unzip the downloaded file:

kallisto index -i 	Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx	Homo_sapiens.GRCh38.cdna.all.fa.gz
```
  
Once the index is created, the transcripts to genes text file must also be compiled. This can be done using a function from **kb_python** called ***create_t2g_from_gtf*** . This requires gtf (gene transfer format) files as input that must be downloaded. 
  
To download gtf files go to ensembl website > human > latest genome assembly > GRCh38 (or latest version) > access the gtf file (Homo_sapiens.GRCh38.110.gtf.gz). Once downloaded, store the gtf file in a specific directory. 

  
```JUPYTER
import kb_python
import re

# For assurance, view the packaged contained within kb_python : 
from kb_python import ref
print(dir(ref))

# This should generate the output: 
['COMBINED_FILENAME', 'FASTA', 'GTF', 'SORTED_FASTA_FILENAME', 'SORTED_GTF_FILENAME', '__builtins__', '__cached__', '__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__', 'concatenate_files', 'create_t2c', 'create_t2g_from_fasta', 'create_t2g_from_gtf', 'download_reference', 'generate_cdna_fasta', 'generate_intron_fasta', 'get_kallisto_binary_path', 'kallisto_index', 'logger', 'logging', 'open_as_text', 'os', 'ref', 'ref_lamanno', 'run_executable', 'sort_fasta', 'sort_gtf', 'tarfile', 'urlretrieve']

# Then to create the transcripts to genes text file :
path_to_gtf = "honours_2023/kallisto_index/Homo_sapiens.GRCh38.110.gtf.gz"
path_to_output = "honours_2023/kallisto_index/transcripts_to_genes"
ref.create_t2g_from_gtf(path_to_gtf, path_to_output, intron=False)

```

