DCDKM.GeneSelection <- function(Mat1, Mat2 = NULL, Cond1type = NULL, Cond2type = NULL
                              , fdr.cut = 0.01, logFC.cut = 1
                              , DEGMethod = "glmLRT"
                              , pipeline="edgeR"
                              , PPI = NULL
                              , PPIcutoff = NULL){


  #### Loading required packages
  if(!require(devtools))
    install.packages("devtools")

  if(!require(AMCBGeneUtils))
    devtools::install_github("AndresMCB/AMCBGeneUtils")

  if(!require(tidyverse))
    install.packages("tidyverse")

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  if(!require(edgeR)){
    BiocManager::install("edgeR")
    library(edgeR)
  }


  if(!is.null(Mat2))
    Features <- intersect(colnames(Mat1),colnames(Mat2)) else
    Features <-colnames(Mat1)
  if(length(Features)==0){
    stop("Error: no common genes between Mat1 and Mat2")
  }

  #### DEG analysis
  if(is.null(Mat2) | is.null(Cond1type) | is.null(Cond2type) ){
    message("DEG analysis is omitted.")
    message("If you want to perform it please be sure to provide MAt2, Cond1type, and Cond2type")
  }  else{
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    if(!require(TCGAbiolinks))
      BiocManager::install("TCGAbiolinks")

    set.seed(1)
    mat1 = t(Mat1)
    # colnames(mat1) <- BRCA_normal$barcode
    mat2 = t(Mat2)
    # colnames(mat2) <- BRCA_PT$barcode

    dataDEGs <- TCGAbiolinks::TCGAanalyze_DEA(mat1 = mat1,
                                              mat2 = mat2,
                                              Cond1type = Cond1type,
                                              Cond2type = Cond2type,
                                              pipeline = pipeline,
                                              metadata = F,
                                              fdr.cut = fdr.cut ,
                                              logFC.cut = logFC.cut,
                                              method = DEGMethod)
    if(any("ID"%in%colnames(dataDEGs)))
      {Features <- dataDEGs$ID}else
        Features <- row.names(dataDEGs)
  }

  #### PPI analysis
  if(is.null(PPIcutoff)){
    message("PPI analysis is omitted.")
    message("If you want to perform it please be sure to provide PPIcutoff >= 1")
  }  else{
    if(is.null(PPI)){
      PPI <- DynamicCancerDriverKM::PPI
      aux1 <- PPI%>%
        dplyr::count(`Input-node Gene Symbol`)
      colnames(aux1) <- c("HGNC.symbol","nIn")
      aux2 <- PPI%>%
        dplyr::count(`Output-node Gene Symbol`)
      colnames(aux2) <- c("HGNC.symbol","nOut")
      PPIdegree <- inner_join(aux1,aux2,by="HGNC.symbol")
      PPIdegree <- PPIdegree%>%
        mutate(nTotal = nIn+nOut)%>%
        arrange(desc(nTotal))

      temp <-
        AMCBGeneUtils::changeGeneId(PPIdegree$HGNC.symbol[PPIdegree$nTotal>PPIcutoff])

      temp <- na.omit(temp$Ensembl.ID)
      Features <- intersect(temp,Features)
    }

  }
  return(Features)
}

