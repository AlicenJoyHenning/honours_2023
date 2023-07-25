# SASCRiP single cell RNA sequencing 
Contents: [Installation](#section-1)

## Installation 
In order for the pipeline to be executed on a device, the following packages must be installed:
+ ***Python*** (greater than version 3.7)
``` python --version ``` if not updated, ``` sudo apt  update ``` / ``` sudo apt install python3 ```
+ ***R*** (greater than version 3.6) ``` R --version ```
+ ***Seurat R package*** R package used for the analysis and visualisation of scRNA-seq data providing a suite of functions including those for quality control and normalisation.
  <pre>
    ```R
    # first check is package is installed:
    library(Seurat)

    # if error message appears, install package: 
    install.packages("Seurat")
    library(Seurat)
    ```
  </pre>
  
+ ***Tidyverse R package*** R package a collection of R packages for data manipulation
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

Once the dependcies are installed and updated, it is now time to install the pipeline itself 
<pre>
  ```Ubuntu command prompt
  # Install package using pip
  pip install SASCRiP 
```
</pre>
