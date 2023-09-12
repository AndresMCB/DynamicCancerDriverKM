
# Summary of experiments
This folder contains a compilation of the results of our experiments.

## ExperimentsSC.rdata
This file contains a list with the results for the 39 explored genes in our experiments. 

The list (named as `results` when you load the data) has 39 elements. Each element corresponds to 1 single experiment (target gene). Name of each element corresponds to the **Ensembl** version of the gene name used as target in the experiment.
  
Each element in the list contains the following objects:
  
  * **topModels**: A list containing the 179 models (out of 47793 explored models) with the best scores towards invariance (_stability score_). A model can be either **_collaborations, interactions, or Main Effect_** model. In these models, terms are the index (starting from 1) of the vector of 179 genes used as model variables in our experiments.

  * **geneScore**: A data frame with 179 rows and 2 columns. The column `features` contains the gene names (Ensemble ID) of the 179 model variables in the order used for modelling (consequently in the order represented in **topModels**). The column `score` contains the _gene score_ that reflects the proportion of stable models that depends on such gene. A large _gene score_ means that such a gene appears in a large proportion of the top models.

  * **InferredDrivers**: A 2 columns dataframe containing genes (`features`) and _gene scores_ (`score`) of the genes with  _gene scores_ grater than 0.
 
  * **formulas**: symbolic representation of the **topModels**.
 
  * **modelScore**: _stability score_ of the top models in the same order as **topModels** and **formulas**. The smaller the _stability score_ the more stable the model towards invariance. 
 
  * **Summary**: A list with 11 elements containing:
    - `CGC.baseline` (element 1):  CGC genes among all 179 genes used as model variables.
    - `CGC_in_top_` (elements 2 to 11): CGC genes among the top 10 to top 80 (elements 2 to 10) (following _gene scores_ in descending order) of **InferredDrivers** and the CGC genes among the **InferredDrivers** (element 11).
    - `summaryTable` (element 12): A dataframe containing the p.value of the hypergeometric test when considering top 10 to top 80 and all **InferredDrivers**.
      
