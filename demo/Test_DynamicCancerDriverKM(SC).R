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

data("GSE75688_TPM", package = "DynamicCancerDriverKM")
data("GSE75688_sample_information", package = "DynamicCancerDriverKM")

SC<- GSE75688_sample_information%>%
  filter(type=="SC")%>%
  select(sample)
GSE75688_TPM_SC <- GSE75688_TPM[SC$sample,,drop=F]
rm(GSE75688_TPM)

#---- remove genes not expressed in at least 20% of the samples ----
GSE75688_TPM_SC <-
  GSE75688_TPM_SC[,colSums(GSE75688_TPM_SC>0)>(0.2*nrow(GSE75688_TPM_SC))
                  , drop=F]

#---- remove version from Ensembl.ID ----
temp <- colnames(GSE75688_TPM_SC)
temp <-  gsub("\\..*","",temp)
colnames(GSE75688_TPM_SC)<- temp
rm(temp)
#---- remove any column with no equivalence un HGNC.symbol ----
library(AMCBGeneUtils)

IdSource <- GeneIdSource(colnames(GSE75688_TPM_SC))
genes <- changeGeneId(colnames(GSE75688_TPM_SC),from = IdSource)
GSE75688_TPM_SC <- GSE75688_TPM_SC[,!is.na(genes$HGNC.symbol),drop=F]

tumour <- GSE75688_sample_information%>%
  filter(type=="SC",index=="Tumor")
non_tumour <- GSE75688_sample_information%>%
  filter(type=="SC",index2=="Stromal")

#mat1 are tumour cells, mat2 are non-tumour
mat1 = GSE75688_TPM_SC[tumour$sample,,drop=F]
mat2 = GSE75688_TPM_SC[non_tumour$sample,,drop=F]

Features <- DCDKM.GeneSelection(Mat1 = mat1, Mat2 = mat2
                              ,Cond1type = "tumour"
                              ,Cond2type = "nonTumour"
                              ,logFC.cut = 1
                              ,PPIcutoff = 1)



# Her2 as covariate for phenopath
binned <- DCDKM.BinTime(Mat1 = rbind(mat1,mat2), covariate = "ENSG00000141736"
                      , Features = Features, nBins = 30)



target <- AMCBGeneUtils::changeGeneId(BRCA.40CD)

target <- intersect(target$Ensembl.ID,colnames(binned$Env1))
results <- vector(mode = "list",length = length(target))
names(results) <- target

counter <- 0
library(tictoc)
for (i in target) {
  tic()
  counter <- counter + 1
  message(paste("Gene ",counter,",", i))
  predictors <- setdiff(Features,i)
  predictors <- intersect(predictors,colnames(binned$Env1))

  # put the target variable "y" as the first feature.
  features <- c(i,predictors)
  # define the number of top models to be considered (k)
  k <- length(features)
  # get the symbolic notation for the models
  testModels <- DCDKM.GetModels(length(features))
  invariantScore<-NULL
  for (j in 1:length(testModels)) {
    tic()
    #Score each model (the lower the score the better towards invariance)
    aux <-  DCDKM.modelScoring(models = testModels[[j]]
                               , binned = binned, parallel = TRUE
                               , features = features
                               , targetIndex = 1, num.folds = 2)
    invariantScore$score <-c(invariantScore$score
                             ,aux$modelScores)
    invariantScore$formulas <-c(invariantScore$formulas
                                ,aux$formulas)
    toc()
  }
  # retrieve the top k models with better "invariant" score
  index <- order(invariantScore$score)
  models <- unlist(testModels, recursive = FALSE)
  topModels <- models[index[1:k]]
  #score genes (the more invariant models the gene appears
  #the higher the score)
  geneScore <- driverScore(topModels, features)
  #Inferred genes as genes with score greater than 0
  InferredDrivers <- geneScore%>%
    filter(score>0)

  #Assemble the results.
  results[[i]]$topModels <- topModels
  results[[i]]$geneScore <- geneScore
  results[[i]]$InferredDrivers <- InferredDrivers
  results[[i]]$formulas <- invariantScore$formulas[index[1:k]]
  results[[i]]$modelScore <- invariantScore$score[index[1:k]]
  results[[i]]$Summary <- overlapCGC(geneScore = geneScore
                                     , top = c(seq(10, 80, 10),nrow(InferredDrivers)))
  toc()
}


####################################################################################
# save the list of inferred driver from Single cell as excel file
# (one target per excel sheet)
temp <-list()
aux <- AMCBGeneUtils::changeGeneId(names(results),from = "Ensembl.ID")[4]

for (i in names(results)) {
  temp[[i]] <-
    cbind(AMCBGeneUtils::changeGeneId(results[[i]]$InferredDrivers[1]
                                      ,from = "Ensembl.ID", to = c("Ensembl.ID","HGNC.symbol"))[2:3]
          ,results[[i]]$InferredDrivers[2])
}

aux <- sapply(aux, function(x){paste("target",x)})
names(temp) <- aux
require(openxlsx)
write.xlsx(temp, file = "supp_Table9 - inferred drivers from single cell.xlsx")

####################################################################################
# save the results of the experiments as a R file for future use

save(results, file = "ExperimentsSC.rdata")



