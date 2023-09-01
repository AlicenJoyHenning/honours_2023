# Summary Visualizations of Gene Ontology Terms With GO-Figure!
Reijnders, M.J. and Waterhouse, R.M., 2021. Summary visualizations of gene ontology terms with GO-Figure!. Frontiers in Bioinformatics, 1, p.6.

_"The Gene Ontology (GO) is a cornerstone of functional genomics research that drives discoveries through knowledge-informed computational analysis of biological data from large-scale assays. Key to this success is how the GO can be used to support hypotheses or conclusions about the biology or evolution of a study system by identifying annotated functions that are overrepresented in subsets of genes of interest. Graphical visualizations of such GO term enrichment results are critical to aid interpretation and avoid biases by presenting researchers with intuitive visual data summaries._

_GO-Figure!, an open-source Python software for producing user-customisable semantic similarity scatterplots of redundancy-reduced GO term lists. The lists are simplified by grouping together terms with similar functions using their quantified information contents and semantic similarities, with user-control over grouping thresholds. Representatives are then selected for plotting in two-dimensional semantic space where similar terms are placed closer to each other on the scatterplot, with an array of user-customisable graphical attributes. GO-Figure! offers a simple solution for command-line plotting of informative summary visualizations of lists of GO terms, designed to support exploratory data analyses and dataset comparisons."_


## [1] Install
First, download this gitlab project or clone ```git clone https://gitlab.com/evogenlab/GO-Figure.git```. The simplest option to then run the GO-Figure! binary files for mac or ubuntu. This is done via the command line using ```./gofigure-ubuntu```. However, in Windows, simply travelling to the directory where the project is saved will work. 
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/d4cb1f95-1e1e-4521-9cc3-16ba0c5b3c7a)


## [2] Load dependencies 
Once GOFigure! project has been cloned, travel to the working directory and run the follwing to ensure it is running properly : 
```python

cd GO-Figure

!python gofigure.py --help

# OUTPUT :
Required arguments:
  -i INPUT, --input INPUT
                        Input file in tabular format with the columns: GO ID +
                        P-value for standard input, GO ID + P-Value +
                        Significant for standard-plus input, TopGO output as
                        an input, or GoStats output as an input. Can use # or
                        ! or % to comment out lines. Change --input_type
                        accordingly. Default input is 'standard'. Example
                        input files are found in the root directory.
  -o OUTPUT, --output OUTPUT
                        Output directory

Optional input file and data handling arguments:
  -a MAX_CLUSTERS, --max_clusters MAX_CLUSTERS
                        Maximum amount of clusters to plot (integer value).
                        Default = 50.
  -j {topgo,gostats,standard,standard-plus}, --input_type {topgo,gostats,standard,standard-plus}
                        Type of input file. Use 'topgo' for standard TopGO
                        input, 'gostats' for standard GOStats input,
                        'standard' for an input file where the first column is
                        the GO ID and the second is the p-value, and
                        'standard-plus' for standard input but with a third
                        column containing a user defined numerical value (for
                        TopGO and GOStats input, this is the 'significant'
                        value). Default = standard
  -n {bpo,mfo,cco,all}, --ontology {bpo,mfo,cco,all}
                        Which ontology to use: biological process ('bpo'),
                        molecular function ('mfo'), cellular component
                        ('cco'), or all ontologies ('all'). Default is all.
  -r REPRESENTATIVES, --representatives REPRESENTATIVES
                        A list of GO terms that have priority for being a
                        representative of a cluster. I.e. if one term in a
                        cluster has priority, that term will always be the
                        representative. If two terms in a cluster have
                        priority, only those two will be considered. Input is
                        a list of GO terms separated by a comma, such as
                        'GO:0000001,GO:0000002'.
  -si SIMILARITY_CUTOFF, --similarity_cutoff SIMILARITY_CUTOFF
                        Similarity cutoff to be used between GO terms, a value
                        between 0 and 1. Default = 0.5.
  -v MAX_PVALUE, --max_pvalue MAX_PVALUE
                        Maximum p-value to consider (floating value). Default
                        = 99999.99999.
  -so {pval,user,user-descending}, --sort_by {pval,user,user-descending}
                        Which values to use for sorting the clusters: 'pval'
                        (p-value) or 'user' (user-defined value) or 'user-
                        descending' (user-defined value descending). Default =
                        pval.
  -su {True,False}, --sum_user {True,False}
                        To sum the user-defined values (column 3) for each
                        member of a cluster. Either 'True' or 'False'. Default
                        = False.
  -to TOP_LEVEL, --top_level TOP_LEVEL
                        Set top level GO terms for clusters as given by the GO
                        DAG (see https://www.ebi.ac.uk/QuickGO). Top level GO
                        terms have to be given separated by comma's, without
                        spaces. E.g.: 'GO:000001,GO:000008'.
  -nc NAME_CHANGES, --name_changes NAME_CHANGES
                        A list of GO terms that will be used as representative
                        of a cluster specifically for naming purposes, but not
                        for internal calculations. This is opposed to the'--
                        representatives option, which provides GO terms to be
                        used as representatives of a cluster both internally
                        and externally. This specific option allows for the
                        renaming of clusters without recalculating the
                        clusters when there is a need to reproduce the
                        original figure. Input is a list of GO terms separated
                        by a comma, such as 'GO:0000001,GO:0000002'.

OPtional :
-p PALETTE, --palette PALETTE
                        Set to any color brewer palette available for
                        MatPlotLib (https://matplotlib.org/3.1.1/tutorials/col
                        ors/colormaps.html). Default = plasma


```
<br>
GOFIgure! works in a Python interface so in a Jupyter Lab notebook the following packages must be loaded before the package can be used for visualisations : 

```python
!pip install numpy
!pip install matplotlib
!pip install seaborn
!pip install scikit-learn
!pip install adjustText
```

## [3] Download updates for GO calculations 
An important aspect of GO-Figure! is the ability to keep up with new versions of data used for its calculations. Updates can be performed by performing the following steps:<br>
+ Download the latest version of the go.obo (<a href="http://geneontology.org/docs/download-ontology/" style="background-color: #6ab5ba; color: white; padding: 10px 20px; border-radius: 5px; text-decoration: none;">GO.obo</a>) and place the downloaded file into the folder from the cloned project named **data**
+ Run the script ```$ python3 scripts/relations.py data/go.obo > data/relations.tsv``` and place the output in the **data** folder 
+  Download the latest version of the GOA UniProt database named ```goa_uniprot_all.gaf.gz``` and place it the **data** folder
+  Run the script ```$ python3 scripts/ics.py data/relations.tab goa_uniprot_all.gaf.gz data/go.obo > data/ic.tsv``` and place the output in the **data** folder  
![image](https://github.com/AlicenJoyHenning/honours_2023/assets/129797527/12846d11-af35-4d67-9889-f583da4409ce)



