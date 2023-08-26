#### ----- Script to test DynamicCancerDriverKM package ------ ####
#
#  This script follows the procedure described in the
#  Briefings in Functional Genomics - Oxford Paper.
#
#
# rm(list = ls())

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("AMCBGeneUtils", quietly = TRUE))
  devtools::install_github(repo = "AndresMCB/AMCBGeneUtils")

if (!requireNamespace("DynamicCancerDriverKM"))
  devtools::install_github(repo = "AndresMCB/DynamicCancerDriverKM")

library(DynamicCancerDriverKM)
library(tidyverse)

# 40 Breast Cancer Drivers used as targets in our experiments
# List from https://www.nature.com/articles/ncomms11479.pdf
BRCA.40CD <-c("SMAD4","USP9X","FOXP1","MEN1", "MLLT4", "TBL1XR1"
              ,"ERBB2","PIK3R1","KDM6A","BRCA1","TP53","ARID1A"
              ,"NF1","SF3B1","MAP2K4","RB1","CDH1","CDKN1B","PIK3CA"
              ,"KMT2C","AGTR2","FOXO3","PTEN","BRCA2","AKT1","TBX3"
              ,"CHEK2","CBFB","NCOR1" ,"KRAS","CDKN2A","RUNX1"
              ,"ZFP36L1","GATA3","MAP3K1","GPS2","CTCF","CTNNA1"
              ,"BAP1","PBRM1")


####---- Load TCGA-BRCA Data -----####
# Dataset downloaded and normalised using TCGAbiolinks (June 2022).
# Samples labelled as Normal Tissue, and samples labelled as Primary Tumour were downloaded.
# Data is included in DynamicCancerDriverKM package as TCGA_BRCA_Normal.rdata and TCGA_BRCA_TP.rdata

data("TCGA_BRCA_Normal", package = "DynamicCancerDriverKM")
data("TCGA_BRCA_PT", package = "DynamicCancerDriverKM")

# Keep only gene expression mat1=normal, mat2 = Primary tumour
mat1 = BRCA_normal[,-c(1:7), drop = F]
row.names(mat1) <- BRCA_normal$barcode
mat2 = BRCA_PT[,-c(1:7), drop = F]
row.names(mat2) <- BRCA_PT$barcode


Features <- DCDKM.GeneSelection(Mat1 = mat1, Mat2 =mat2, Cond1type = "normal",Cond2type = "PT"
                              , PPIcutoff = 4)

# Her2 as covariate for phenopath
binned <- DCDKM.BinTime(Mat1 = rbind(mat1,mat2), covariate = "ENSG00000141736"
                      , Features = Features)

# Run experiments only for those targets that are present in the dataset.
target <- AMCBGeneUtils::changeGeneId(BRCA.40CD)
target <- intersect(target$Ensembl.ID,colnames(binned$Env1))

results <- vector(mode = "list",length = length(target))

#i <- target[1]
#library(tictoc)
for (i in target) {
  #tic()
  predictors <- setdiff(Features,i)
  predictors <- intersect(predictors,colnames(binned$Env1))
  features <- c(i,predictors)
  k <- length(features)

  testModels <- DCDKM.GetModels(length(features), option = 3)
  models <- c(testModels$models.1, testModels$models.2, testModels$models.3)

  invariantScore<-NULL
  for (j in 1:3) {
    #tic()
    aux <-  DCDKM.modelScoring(models = testModels[[j]]
                         , binned = binned, parallel = TRUE
                         , features = features
                         , targetIndex = 1, num.folds = 2)
    invariantScore$score <-c(invariantScore$score
                             ,aux$modelScores)
    invariantScore$formulas <-c(invariantScore$formulas
                                ,aux$formulas)
   # toc()
  }


  index <- order(invariantScore$score)
  topModels <- models[index[1:k]]

  geneScore <- driverScore(topModels, features)

  DCDK <- geneScore%>%
    filter(score>0)
  #performance.CGC(geneScore = geneScore, top = seq(50, 285, 5))

  #formulas <- sapply(topModels,getFormula, features)
  #topModels$InvariantScore <- invariantScore$score[index[1:k]]

  results[[i]]$topModels <- topModels
  results[[i]]$geneScore <- geneScore
  #  results[[i]]$DCDk <- DCDK
  results[[i]]$formulas <- invariantScore$formulas[index[1:k]]
  #  results[[i]]$fits <- invariantScore$fits[index[1:k]]
  results[[i]]$InvariantScore <- invariantScore$score[index[1:k]]
  results[[i]]$Summary <- performance.CGC(geneScore = geneScore
                                          , top = seq(50, 250, 50))
  #toc()
}





