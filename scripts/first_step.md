```python
import SASCRiP
```


```python
from SASCRiP import sascrip_functions
```


```python
sascrip_functions.install_R_packages()
```

    Error in contrib.url(repos, type) : 
      trying to use CRAN without setting a mirror
    Calls: install.packages -> startsWith -> contrib.url
    Execution halted
    



```python
sascrip_functions.install_R_packages()
```

    Error in contrib.url(repos, type) : 
      trying to use CRAN without setting a mirror
    Calls: install.packages -> startsWith -> contrib.url
    Execution halted
    



```python
import os 
```


```python
os.listdir("Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/")
```




    ['bamtofastq_S1_L007_I1_003.fastq.gz',
     'bamtofastq_S1_L007_R2_002.fastq.gz',
     'bamtofastq_S1_L007_I1_007.fastq.gz',
     'bamtofastq_S1_L007_I1_001.fastq.gz',
     'bamtofastq_S1_L007_R2_006.fastq.gz',
     'bamtofastq_S1_L007_R1_005.fastq.gz',
     'bamtofastq_S1_L007_I1_002.fastq.gz',
     'bamtofastq_S1_L007_R1_007.fastq.gz',
     'bamtofastq_S1_L007_R1_006.fastq.gz',
     'bamtofastq_S1_L007_I1_006.fastq.gz',
     'bamtofastq_S1_L007_I1_004.fastq.gz',
     'bamtofastq_S1_L007_R1_003.fastq.gz',
     'bamtofastq_S1_L007_R2_005.fastq.gz',
     'bamtofastq_S1_L007_R1_001.fastq.gz',
     'bamtofastq_S1_L007_R2_007.fastq.gz',
     'bamtofastq_S1_L007_R1_002.fastq.gz',
     'bamtofastq_S1_L007_R1_004.fastq.gz',
     'bamtofastq_S1_L007_R2_001.fastq.gz',
     'bamtofastq_S1_L007_R2_003.fastq.gz',
     'bamtofastq_S1_L007_I1_005.fastq.gz',
     'bamtofastq_S1_L007_R2_004.fastq.gz']




```python
sascrip_functions.kallisto_bustools_count(
    "Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/"
```


```python
# set parameters 
list_of_fastqs = [
    "Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/bamtofastq_S1_L007_R1_001.fastq.gz", 
    "Alicen/raw_files/my_fastq_files/human_blood_ifnalpha/bamtofastq_S1_L007_R2_001.fastq.gz"
]
single_cell_technology = "10Xv3"
output_directory_path = "Alicen/test_files"
species_index = "Alicen/kallisto_index.idx"
species_t2g = "Alicen/transcripts_to_genes.txt"

```


```python
sascrip_functions.kallisto_bustools_count(
    list_of_fastqs = list_of_fastqs,
    single_cell_technology = single_cell_technology,
    output_directory_path = output_directory_path, 
    species_index = species_index,
    species_t2g = species_t2g
)



                                          
```




    {'transcript_to_genes': 'Alicen/transcripts_to_genes.txt',
     'ec_mapping_file': 'Alicen/test_files/Count_analysis/matrix.ec'}




```python

```
