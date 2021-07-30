# Barcode Fusion Genetics related codes

This page explains how to install and execute the BFG-PCA codes.
Please go to the [BFG-PCA wiki (still under construction)](https://github.com/DanYamamotoEvans/BFG-PCA/wiki) for more information regarding each step of analysis.
If you have any questions, pleaes post a question in the Discussions. 

Please make sure you have appropriate Python and pip before starting.
```sh
Python version >=3.5
pip    version >= 1.1.0
```

Dependencies :
```sh
numpy  version >=1.19 
tqdm   version >=2.2.4
```
To install these pakcages, first clone this repository by
```sh
git clone https://github.com/DanYamamotoEvans/BFG-PCA.git
```

Next, go to the location of the BFG-PCA folder in the terminal, and install the dependencies by
```sh
pip install .
```

Other core programs to install:
- [Jupyter-notebook](https://jupyter.org/install)
```sh
pip install jupyterlab
```
- Commandline BLAST+


Follow the [instruction manual](https://www.ncbi.nlm.nih.gov/books/NBK569861/) for installation.
Set the PATH of the binary file.

Execute the following to see if installation is complete.
```sh
blastn -help
```
    
## Overview
This script was built to perform experimental plans and data analysis for Barcode Fusion Genetics screenings. There are four main steps in this suit, which I have prepared jupyter-notebooks for each.

- Monte-Carlo simulation
- Barcode calling
- Normalization
- Performance measure
- (Visualization, You will need to install [R](https://cran.r-project.org/).)

### Monte-carlo simulation of BFG screening proccess
Since BFG screenings have multiple sampling steps while handling a complex pool of strains, we suimulate the sampling process with a Monte-Carlo simulation. This notebook follows the procedures of BFG screenings, and allows the user to estimate the nessesary paramaters for sampling. 

### Barcode calling
This notebook creates the BLAST databse based on your barcode database file, and performs BLAST on the fastq files you provide. The results will be parsed to combine the count data. Please look at the wiki for input file specifications.

### Normalization
This notebook normalizes the count data, and compute raw PPI signals based on the barcode counts in the control condition and auto-activity level.
It will also output some csv files for plotting stats.

### Performance measure
Based on the normalized scores, this notebook computes the agreement against the BioGRID database. It willcomute the agreement for various percentile values of the PPI scores generated from multiple replicates.


### Visualization
This notebook will help you plot the data you obtained from the BFG screening. Please go to the visualization page in the wiki for more details on the plots.

### References
- [Yachie _et al_, 2016](https://www.embopress.org/doi/full/10.15252/msb.20156660) / Initial report of BFG. The codes here were built based on perl scripts provided from [Dr. Nozomu Yachie](http://yachie-lab.org/?nozomuyachie).
- [Evans-Yamamamto _et al_, 2021 (Preprint)](https://www.biorxiv.org/content/10.1101/2021.07.27.453987v1) / This repositry was created in part of this work to make BFG-PCA analysis accessible. 
