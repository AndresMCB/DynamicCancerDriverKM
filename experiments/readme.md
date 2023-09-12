
# Summary of experiments
This folder contains a compilation of the results of our experiments.

## ExperimentsSC.rdata
This file contains a list with the results for the 39 explored genes in our experiments. 
* The list (named as `results` when you load the data) has 39 elements. Each element corresponds to 1 single experiment (target gene). Name of each element corresponds to the **Ensembl** version of the gene name used as target in the experiment.
* Each element in the list contains the following objects:
  - **topModels**: A list containing the 179 models (out of 47793 explored models) with the best scores towards invariance. A model can be either **_collaborations, interactions, or Main Effect_** model. 
