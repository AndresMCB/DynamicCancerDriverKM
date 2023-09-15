
# Summary of experiments
This folder contains a compilation of the results of our experiments. In both cases, loading the data will create an R list named `results`. Elements in `results` are named as the Ensembl ID of the target gene used for that especific experiment. If you want to rename the elements to its HGNC symbol (gene symbol) you can run the following code after loading `results`:

```R
# change names of the targets to "target" + HGNC symbol for the sake of simplicity
aux <- AMCBGeneUtils::changeGeneId(names(results),from = "Ensembl.ID")[4]
aux <- sapply(aux, function(x){paste("target",x)})
names(temp) <- aux
```
## ExperimentsBulk.rdata

This file contains a list (`results`) with 39 elements corresponding to the results for the 39 explored genes in our experiments when using the TCGA-BRCA dataset (TCGA files in [data](\data) folder).  Name of each element of the `results` list corresponds to the **Ensembl** ID of the gene used as target in that specific experiment.

## ExperimentsSC.rdata

This file contains a list (`results`) with 39 elements corresponding to the results for the 39 explored genes in our experiments when using the single cell RNA sequencing data from NCBI GEO database, accession [GSE75688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688).  Name of each element corresponds to the **Ensembl** ID of the gene used as target in the experiment.

##  Structure of the `results` lists

Each element in the `results` list (both bulk data and single cell data) contains the following objects:
  
  * **topModels**: A list containing the models (**TCGA-BRCA**: 733 out of 804834 explored models. **GSE75688**: 179 out of 47793 explored models) with the best scores towards invariance (_stability score_). A model can be either **_collaborations, interactions, or Main Effect_** model. In these models, terms are the index (starting from 1) of the vector of genes used as model variables in our experiments.

  * **geneScore**: A data frame with 2 columns. The column `features` contains the gene names (Ensemble ID) of all model variables in the order used for modelling (consequently, in the order represented in **topModels**). The column `score` contains the _gene score_ that reflects the proportion of stable models that depends on such gene. A large _gene score_ means that such a gene appears in a large proportion of the top models.

  * **InferredDrivers**: A 2 columns dataframe containing genes (`features`) and _gene scores_ (`score`) of the genes with  _gene scores_ greater than 0.
 
  * **formulas**: A character array with the symbolic representation of the **topModels**.
 
  * **modelScore**: _stability score_ of the top models in the same order as **topModels** and **formulas**. The smaller the _stability score_ the more stable the model towards invariance. 
 
  * **Summary**: A list with 8 elements (**TCGA-BRCA**) or 11 elements (**GSE75688**) containing:
    - `CGC.baseline` (element 1):  CGC genes among all genes used as model variables.
      
    - `CGC_in_top_` (**TCGA-BRCA**: elements 2-7, **GSE75688**: elements 2-11): CGC genes among the different explored tops of **InferredDrivers** (following _gene scores_ in descending order). The last element (**TCGA-BRCA**: element 7, **GSE75688**: element 11) corresponds to the CGC genes among all **InferredDrivers**.

    - `summaryTable` (**TCGA-BRCA**: element 8, **GSE75688**: element 12): A dataframe containing the p.value of the hypergeometric test when considering all explored tops of **InferredDrivers**.
      
