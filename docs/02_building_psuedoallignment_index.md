
# Building the reference transcriptome input for kallisto 

To run the pseudo alignment tool (kallisto), an index of the reference transcriptome is needed. Although the SASCRiP function **kallisto_bustools** is able to do this automatically by changing some parameters, I needed to know how to do this manually.<br> 
## [1] Install Dependencies<br>
For this process, python needs to be installed along with the package manager, pip. Download the latest version of python (<a href="https://www.python.org/downloads/" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">PYTHON</a>) and download the file here: (<a href="https://bootstrap.pypa.io/get-pip.py" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">PIP</a>). Ensure Python is installed from the cmd using ```python --version``` before installing pip ```python get-pip.py``` which can be verified afterwards as well using ```pip --version```. Next, JupyterLab must be installed which can now be done from the cmd ```pip install jupyterlab``` and opened  ```jupyter lab```. Note that this takes you to a Google Chrome page with Jupyter ready to be used, this is where the second portion of this tutorial needs to be completed.

This process also requires Kallisto to be installed (<a href="https://github.com/pachterlab/kallisto/releases" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">KALLISTO</a>). From the options, choose and install the one compatible with your device and to use it, travel to the directory where it is saved ```cd kallisto_windows-v0-50.0\kallisto``` : 
+ kallisto_linux-v0.50.0.tar.gz
+ kallisto_mac-v0.50.0.tar.gz
+ kallisto_mac_m1-v0.50.0.tar.gz
+ kallisto_windows-v0.50.0.zip
<br>
From the command prompt, 

## [2] Download the necessary 
Once Kallisto is installed and ready to be used , the full transcriptome from Ensembl (files ending in cdna.all.fa.gz) must be downloaded. To build the human transcriptome index, first download the transcriptome, which is available under cDNA on the Ensembl website (<a href="http://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">ENSEMBL</a>) and execute the following in the command prompt : 

  
```command promt 
# Download the full transcriptome from ensemble : 
curl -O ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

<div style="text-align:center">
    <img src="https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/7c5f8b9b-e275-4dd4-b79f-1ddc0e55b37f" alt="image">
</div>
_NOTE: you could, of course, do this by clicking the download button when you travel to the website_ <br><br><br>
After downloading the full transcriptome file, you now need to build the index :<br>

```
# Run kallisto index. kallisto will work on .fa and .fz.gz files so there is no need to unzip the downloaded file:

kallisto index -i 	Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx	Homo_sapiens.GRCh38.cdna.all.fa.gz
# -i Homo_sapiens.GRCh38.cdna.all.release-94_k31.idx: specifies the output index file name
# Homo_sapiens.GRCh38.cdna.all.fa.gz: input transcriptome reference file where -i flag indicates that you want to build an index from this file
# OUTPUT :
[build] building MPHF
[build] creating equivalence classes ...
[build] target de Bruijn graph has k-mer length 31 and minimizer length 23
[build] target de Bruijn graph has 936536 contigs and contains 108619921 k-mers 

```
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/ecef3b27-fc09-4847-801e-42bc643877d6)<br><br>




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

ALternatively using R (<a href="https://bioconductor.org/packages/devel/bioc/vignettes/BUSpaRse/inst/doc/tr2g.html" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">BUSpaRse</a>) 
