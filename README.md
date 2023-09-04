# DynamicCancerDriverKM
DynamicCancerDriverKM package contains functions to identify genes interacting and collaborating to drive alterned core biological processes during cancer progression. Our package test causality in the setup of dynamical system models. Formally, our method test the causal structure of gene collaborations and interactions durng cancer development to identify cancer drivers and collaborative gene networks. 

### Note
Please find the results of our experiments in the folder [inst/Supplementary/](inst/Supplementary/). You can also access them from R by using the name of the file. For example:

```R
aux <- system.file("Supplementary/"
                   , "supplementary table 1 hypertest for 39 target genes(Bulk).csv"
                   , package = "DynamicCancerDriverKM")
Summary_bulk <- read.csv(aux)
```

## Introduction 
Our method takes gene expression data from cross-sectional studies. 

We use the pseudotime (either provided or inferred by using PhenoPath) and the covariate provided to
find a critical turning point in the trajectory along the pseudotime. We
name this critical point as the "event".

We applied our dynamic cancer drivers approach to a single cell RNA
sequencing dataset (NCBI GEO database, accession [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)) ([Chung
et al., 2017](https://www.nature.com/articles/ncomms15081)), and the cancer genome atlas breast cancer dataset ([TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)).
Experiments implemented in our paper can be found as follows:
1. [demo/Test_DynamicCancerDriver(SC).R](demo/Test_DynamicCancerDriver(SC).R): Drivers inferred from a pre-processed Single Cell data, (GSE75688)
2. [demo/Test_DynamicCancerDriver(Bulk).R](demo/Test_DynamicCancerDriver(Bulk).R): Drivers inferred from the TCGA-BRCA dataset.

## Installation 
DynamicCancerDriver runs in the R statistical computing environment.

R (>=4.1.0), devtools(>=2.4.2), Rtools (>=4.0), Bioconductor (>=3.14), phenopath (tested with v. 1.18.0), CausalImpact(>= 1.2.7), and
 tidyverse(>= 1.3.1) are  required.
We also use some utilities from another of our packages ([AMCBGeneUtils](https://github.com/AndresMCB/AMCBGeneUtils)).

1. Please download and install Rtools 4.0 from https://cran.r-project.org/bin/windows/Rtools/, remove the incompatible version from your PATH.

2. Please install devtools 

```R
install.packages("devtools")
```

3. Install DynamicCancerDriver package from github repository 
```R
devtools::install_github('AndresMCB/DynamicCancerDriver')
```

4. Please install DynamicCancerDriver additional packages required as follows: 
```R
devtools::install_github('AndresMCB/AMCBGeneUtils')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phenopath")

```
## Documentation 
Detailed information about the functions implemented in PTC can be found in the [user manual](DynamicCancerDriver_1.4.1.pdf)

Please find the datasets employed in our paper in the folder [data](data/)
