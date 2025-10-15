# DynamicCancerDriverKM

**DynamicCancerDriverKM** identifies cooperative cancer driver genes using *dynamic causal inference*.   It integrates **pseudotime analysis** to order samples along a temporal trajectory and models gene interactions through **causal kinetic models**, an extension of structural causal models to dynamic systems.  

By combining **biological knowledge** with **temporal ordering** (causes preceding effects), the framework expands causal discovery beyond static data, enabling the identification of both **mutated** and **non-mutated** driver genes from **bulk** and **single-cell** datasets.

## Introduction

Our method takes gene expression data from cross-sectional studies and integrates a temporal dimension along cancer progression using pseudotime analysis. This allows testing the **causality of candidate genes** within a dynamic framework, rather than relying on static data. We apply **causal kinetic models** to identify gene interactions and collaborations that drive cancer development.

We have evaluated our method on both **single-cell** and **bulk sequencing** datasets of breast cancer, including the [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688) dataset ([Chung *et al.*, 2017](https://www.nature.com/articles/ncomms15081)) and the [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) cohort. The results show that our approach systematically identifies bona fide driver genes and detects sets of genes strongly linked to cancer progression (including **mutated and non-mutated** drivers) providing a more comprehensive view of cancer development.

As part of this work, we have also developed an **in silico gene expression simulator**, which generates synthetic dynamic datasets to evaluate model performance and assess robustness under controlled conditions. This addition was inspired by a valuable suggestion from our reviewers, and it significantly enhances our capacity to test and benchmark causal inference strategies under diverse simulated biological scenarios.

## Experiments

Scripts with the experiments implemented in our paper can be found in the folder [demo](\d) as follows:
1. [demo/Test_DynamicCancerDriverKM(SC).R](demo/Test_DynamicCancerDriverKM(SC).R): Experiments using a pre-processed Single Cell data, (GSE75688)
2. [demo/Test_DynamicCancerDriverKM(Bulk).R](demo/Test_DynamicCancerDriverKM(Bulk).R): Experiments using the TCGA-BRCA dataset.
3. [Simulations](Simulations): Folder with our in-Silico experiments.In includes the scripts to generate syntethic data, and evaluation of our method in such data. It also includes the results reported in our paper.

> [!IMPORTANT]  
> Running the demo scripts can be time consuming (especially for the TCGA dataset). For personal computers we recommend to test our method by using one target gene at the time. Alternatively, you can find the outcomes of our demo scripts in the (results of our experiments as R data lists) in the [experiments](experiments/) folder.
> 

## Installation 
DynamicCancerDriverKM runs in the R statistical computing environment.

R (>=4.1.0), devtools(>=2.4.2), Rtools (>=4.0), Bioconductor (>=3.14), phenopath (tested with v. 1.18.0), and
 tidyverse(>= 1.3.1) are  required.
We also use some utilities from another of our packages ([AMCBGeneUtils](https://github.com/AndresMCB/AMCBGeneUtils)).

1. Windows only: please download and install Rtools from https://cran.r-project.org/bin/windows/Rtools/, remove any incompatible version from your PATH.

1. Please install devtools

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

* Supplementary tables described in our paper (including the full list of inferred drivers of each experiment) can be found in the folder [supplementary](supplementary/)
