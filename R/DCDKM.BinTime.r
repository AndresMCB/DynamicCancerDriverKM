DCDKM.BinTime <- function(Mat1, Mat2 = NULL, covariate, Features
                        , elbo_tol = 1e-03, nBins=50, pTime = NULL){

  binningGE <- function(GE,nBins, pTime){
    GE <- GE[order(pTime),,drop=F]
    pTime <- pTime[order(pTime),,drop=F]
    indexes <- matrix(nrow=nBins,ncol=2)
    indexes[,1] <- as.integer(seq(1,nrow(GE)+1,length.out=nBins+1))[-(nBins+1)]
    indexes[,2] <- as.integer(seq(1,nrow(GE)+1,length.out=nBins+1)-1)[-1]

    binned <- apply(indexes,MARGIN = 1
                    ,function(index,c)
                    {
                      colMeans(c[index[1]:index[2],])
                    }
                    ,cbind(pTime,GE))
    binned <- t(binned)
    return(list(binnedGE = binned[,-1,drop=F], binTime = binned[,1,drop=F] ))
  }

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  if (!require("phenopath",  quietly = TRUE)){
    BiocManager::install("phenopath")
    library("phenopath")
  }

  if (!require("edgeR")){
    BiocManager::install("edgeR")
    library("edgeR")
  }

  # if gene expression of 2 conditions is provided,
  # bind them in a single gene expression matrix G
  if(!is.null(Mat2)){
    G <- rbind(Mat1,Mat2)
  } else
  {
    G <- Mat1
  }

  findPtime <- TRUE
  if(!is.null(pTime)){
    pTime <- as.matrix(temp, ncol=1)
    if(nrow(pTime) != nrow(G)){
      message(paste("pTime provided is not suitable."
      ,"pTime should be a vector with length equal to"
      ,"number of rows (samples)")
      )
      message("using Phenopath to calculate a pseudotime order")
    }else{
      findPtime <- FALSE
      row.names(pTime) <- row.names(exprs_obj)
    }
  }
  if(findPtime){
    covariate <-intersect(covariate, colnames(G))
    if(length(covariate)==0)
      stop("Error: Covariate was not found in G")

    # Find the number of bins
    if(!is.numeric(nBins) || nBins<10)
      stop("Error: number of Bins (nBins) needs to be a
         positive number greater than 10")

    # Find the number of samples per bin
    spBin <- nrow(G)/nBins
    if((spBin)<2)
      stop("number of bins (nBins) is too large.
         At least 2 samples per bin are required.")
      exprs_obj <-data.matrix(G)
      Features <- intersect(Features, colnames(exprs_obj))
      aux <- setdiff(Features, covariate)
      exprs_obj <- log2(exprs_obj[,aux,drop=F]+1)

      PhenoP <- phenopath(exprs_obj = exprs_obj
                          , x = data.matrix(G[,covariate, drop=F])
                          , elbo_tol = elbo_tol)
      plot_elbo(PhenoP)
      pTime <- as.matrix(trajectory(PhenoP),ncol=1)
      row.names(pTime) <- row.names(exprs_obj)
  }

    #ascending order
    pOrder <- row.names(pTime[order(pTime), ,drop=F])
    pTime1 <- pTime[pOrder[seq(1,length(pOrder),2)],,drop=F]
    pTime2 <- pTime[pOrder[seq(2,length(pOrder),2)],,drop=F]

    Env2 <- G[row.names(pTime2),, drop = F]
    Env1 <- G[row.names(pTime1),, drop = F]

  # binning Gene Expression

  Env1 <- binningGE(GE=Env1,nBins, pTime = pTime1)
  Env2 <- binningGE(GE=Env2,nBins, pTime = pTime2)

  binTime1 <- Env1$binTime
  binTime2 <- Env2$binTime
  binTime <-rowMeans(cbind(binTime1,binTime2))
  return(list(Env1 = Env1$binnedGE,Env2=Env2$binnedGE, binTime = binTime
              , pTime = pTime))


}
