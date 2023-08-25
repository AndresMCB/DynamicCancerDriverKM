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


