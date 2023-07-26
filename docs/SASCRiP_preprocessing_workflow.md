# SASCRiP single cell RNA sequencing 
Contents: 
[(i) Installing dependencies](#section-1) 
[(ii) Installing SASCRiP](#section-2) 
[(iii) Preparing fastq files for SASCRiP](#section-3)



## (i) Installing dependencies 
In order for the pipeline to be executed on a device, the following packages must be installed:
+ ***Python*** (greater than version 3.7)
``` python --version ``` (in my case it was Python 3.10.11) if not updated, ``` sudo apt  update ``` / ``` sudo apt install python3 ```
+ ***R*** (greater than version 3.6) ``` R --version ``` (in my case it was R version 4.2.3)
+ ***Seurat R package*** (in my case it wasn't installed yet) R package used for the analysis and visualisation of scRNA-seq data providing a suite of functions including those for quality control and normalisation.
  <pre>
    ```R
    # first check is package is installed:
    library(Seurat)

    # if error message appears, install package: 
    install.packages("Seurat")
    library(Seurat)
    ```
  </pre>
  
+ ***Tidyverse R package*** (in my case this was installed already) R package a collection of R packages for data manipulation
 <pre>
    ```R
    library(tidyverse)
   
   # Otherwise:
   
    install.packages("tidyverse")
    library(tidyverse)
    ```
  </pre>
  
+ ***kb-python*** (**kallisto|bustools**) for pseudo-alignment and quantification capabilities of Kallisto and the single-cell processing utilities of bustools
<pre>
  ```Ubuntu command prompt
  # Check the version available 
  kb --version 

  # If not installed yet, first update **pip** 
# pip is a package manager for python (Pip Installs Packages)
  
  python3 -m pip install --upgrade pip
  
# Install **kb-python** using **pip** 
  sudo python3 -m pip install kb-python

  # To confirm the installation check the version available 
  kb --version 
```
</pre>


## (ii) Installing SASCRiP  
Once the dependcies are installed and updated, it is now time to install the pipeline itself 
<pre>
  ```Ubuntu command prompt
  # Install package using pip
  pip install SASCRiP 

  # Install additional R packages required for SASCRiP using the function 'install_R_packages'
  import SASCRiP 
  from SASCRiP import sascrip_functions 

  sascrip_functions.install_R_packages()
  
```
</pre>

## (iii) Preparing fastq files for SASCRiP 


