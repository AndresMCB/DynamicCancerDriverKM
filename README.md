# DynamicCancerDriverKM
DynamicCancerDriverKM package contains functions to identify genes interacting and collaborating to drive altered core biological processes during cancer progression. Our package test causality in the setup of dynamical system models. Formally, our method test the causal structure of gene collaborations and interactions during cancer development to identify cancer drivers and collaborative gene networks. 

## Introduction 
Our method takes gene expression data from cross-sectional studies. The method integrates the temporal dimension of the data along the cancer progression and provides a way to test for causality of candidate genes on cancer. We have applied our method to single-cell and bulk sequencing datasets of breast cancer. The evaluation results show that our method systematically identifies bona fide driver genes and detects sets of genes strongly linked to cancer progression. The results suggest that
our method can discover mutated and non mutated drivers of cancer to provide a comprehensive view of cancer development.

We applied our dynamic cancer drivers approach to a single cell RNA sequencing dataset (NCBI GEO database, accession [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688)) ([Chung et al., 2017](https://www.nature.com/articles/ncomms15081)), and the cancer genome atlas breast cancer dataset ([TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)).

## Experiments
Scripts with the experiments implemented in our paper can be found as follows:
1. [demo/Test_DynamicCancerDriverKM(SC).R](demo/Test_DynamicCancerDriverKM(SC).R): Experiments using a pre-processed Single Cell data, (GSE75688)
2. [demo/Test_DynamicCancerDriverKM(Bulk).R](demo/Test_DynamicCancerDriverKM(Bulk).R): Experiments using the TCGA-BRCA dataset.

Results of our experiments (as R data lists) can be found in the folder [experiments](experiments/)
1.

## Installation 
DynamicCancerDriver runs in the R statistical computing environment.

R (>=4.1.0), devtools(>=2.4.2), Rtools (>=4.0), Bioconductor (>=3.14), phenopath (tested with v. 1.18.0), and
 tidyverse(>= 1.3.1) are  required.
We also use some utilities from another of our packages ([AMCBGeneUtils](https://github.com/AndresMCB/AMCBGeneUtils)).

1. Please download and install Rtools 4.0 from https://cran.r-project.org/bin/windows/Rtools/, remove the incompatible version from your PATH.

2. Please install devtools 

```R
install.packages("devtools")
```

3. Install DynamicCancerDriverKM package from github repository 
```R
devtools::install_github('AndresMCB/DynamicCancerDriverKM')
```

4. Please install additional packages required for DynamicCancerDriverKM as follows: 
```R
devtools::install_github('AndresMCB/AMCBGeneUtils')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR","phenopath", "TCGAbiolinks"))

```
## Documentation 
* Please find the datasets employed in our paper in the folder [data](data/)

* Supplementary tables described in our paper can be found in the folder [supplementary](supplementary/)
