summariseResults <- function(resultList, top = seq(50, 150, 50) ){
  qwer <- lapply(resultList
                 , function(GeneResult){
                   topModels <- GeneResult$topModels
                   features <- GeneResult$geneScore[["features"]]
                   geneScore <- driverScore(topModels, features)
                   outcome <- performance.CGC(geneScore = geneScore
                                              , top = top)
                   outcome$geneScore <- geneScore
                   return(outcome)
                 })
  newNcol <- 2*ncol(qwer[[1]][["summaryTable"]])+4
  allResults2 <- matrix(nrow = length(qwer)
                        , ncol = newNcol)
  for (i in 1:length(qwer)) {
    allResults2[i,1] <- names(qwer[i])
    allResults2[i,2] <- length(qwer[[i]][[1]])
    allResults2[i,3] <- nrow(qwer[[i]][["geneScore"]])
    allResults2[i,4] <- sum(qwer[[i]][["geneScore"]]["score"]>0)
    aux <- performance.CGC(geneScore = qwer[[i]][["geneScore"]]
                           , top = as.numeric(allResults2[i,4]))
    allResults2[i,5:6] <- data.matrix(cbind(length(aux[[2]]),aux[[3]][2]))
    for (j in 1:(ncol(qwer[[i]]$summaryTable)-1)) {
      allResults2[i,2*j+5] <- length(qwer[[i]][[(j+1)]])
      allResults2[i,2*j+6] <-
        data.matrix(qwer[[i]]$summaryTable[(j+1)])
    }
  }
  colnames(allResults2) <- paste0("c",1:ncol(allResults2))
  colnames(allResults2)[1:6] <- c("target","CGC_baseline", "nVars", "nDCDK", "top_nDCDK"
                                  ,"pval_nDCDK")
  colnames(allResults2)[seq(7,newNcol,2)] <- paste0("top_",top)
  colnames(allResults2)[seq(8,newNcol,2)] <- paste0("pval_top_",top)
  allResults2 <- cbind(AMCBGeneUtils::changeGeneId(allResults2[,1]
                                                   , from = "Ensembl.ID"
                                                   ,to = "HGNC.symbol")[2]
                       , allResults2)
  #View(allResults2)
  return(allResults2)
}



getFormula <- function(model, features, targetIndex){
  f <- paste(paste0("d_",features[targetIndex]),"~ ")
  for (i in 1:length(model)) {
    if(i>1)
      f <- paste(f,"+")
    f <- paste(f,paste(features[model[[i]]],collapse = ":"))
  }
  f <- paste(f,"- 1")
  return(f)
}


performance.CGC<-function(geneScore = NULL
                          , top = c(50,100,150,200,250)
){
  aux <- list()
  features <-geneScore$features
  aux$CGC.baseline = intersect(geneScore$features ,CGC.driverNames$Ensembl.ID)
  geneScore <- geneScore %>%
    arrange(desc(score))%>%
    filter(score>0)


  summaryTable <- data.frame(CGC.baseline = length(aux$CGC.baseline))

  for (i in top) {
    varName <- paste0("top_", i)
    aux[[varName]] <- intersect(geneScore$features[1:i]
                                ,CGC.driverNames$Ensembl.ID)

    p.val <- 1 - phyper(q = length(aux[[varName]])-1
                        , m = length(aux$CGC.baseline)
                        , n = length(features)-length(aux$CGC.baseline)
                        , k = i, lower.tail = T, log.p = FALSE)
    varName <- paste0("p.val_top_", i)
    summaryTable[[varName]] <- p.val
  }
  aux$summaryTable <- summaryTable
  return(aux)
}

##############
##############




DCDKM.modelScoring <- function(models = NULL, binned, features, targetIndex, parallel = T
                         , num.folds = 2, score.type = "mean_absolute"){

  if(!require(fda, quietly = T)){
    install.packages("fda")
    library(fda)
  }

  if(!require(cvTools, quietly = T)){
    install.packages("cvTools")
    library(cvTools)
  }

  if(!require(quadprog, quietly = T)){
    install.packages("quadprog")
    library(quadprog)
  }


  k <- length(features)
  num.folds <- 2
  integrated.model <- T
  num.models <- length(models)

  Y <- cbind(binned$Env1[,features[targetIndex],drop=F]
             , binned$Env2[,features[targetIndex],drop=F])
  L <- nrow(Y)

  dYlist <- vector("list", 2)
  RSS_A <- vector("numeric", 2)
  lambda <- vector("numeric", 2)

  scores <-  vector("numeric", num.models)
  ## Find data fit A "Data Fit A calculates a smoothing spline to the data
  ## using all realization from the same experiment"
  for (i in 1:ncol(Y)) {
    #fit <- CausalKinetiX::
    fit <- constrained.smoothspline(Y[,i],
                                    binned$binTime,
                                    pen.degree=2,
                                    constraint="none",
                                    times.new=0,
                                    num.folds=num.folds,
                                    lambda="optim")
    dYlist[[i]] <- fit$smooth.deriv
    # compute differences for integrated model fit
    if(integrated.model){
      dYlist[[i]] <- diff(Y[,i,drop=FALSE])
      colnames(dYlist[[i]]) <- paste0("d_",colnames(dYlist[[i]]))
    }
    # fit the reference model and compute penalty par and RSS
    lambda[i] <- fit$pen.par
    RSS_A[i] <- sum(fit$residuals^2)
  }

  if(parallel){
    if(!require(parallel))
      install.packages("parallel")

    copies_of_r <- detectCores(logical = T)-1
    cl <- makeCluster(copies_of_r)
    include.intercept = F
    clusterExport(cl, c("models","binned", "features", "integrated.model"
                        , "L", "dYlist", "num.folds", "Y", "lambda"
                        , "score.type", "getFormula", "constrained.smoothspline"
                        )
                  ,envir = environment())
    modelScores <-
      parSapply(cl,1:num.models, simplify = F#T
                , function(modelIndex){
                  RSS_B <- vector("numeric", 2)
                  Xlist <- vector("list", 2)
                  # fitLM <- vector("list", 2)
                  # fitSpline <- vector("list", 2)
                  ## collect predictors X
                  index <- unique(unlist(models[[modelIndex]]))
                  Xlist[[1]] <- binned$Env1[,features[index],drop=F]
                  Xlist[[2]] <- binned$Env2[,features[index],drop=F]

                  formula <- getFormula(model = models[[modelIndex]]
                                        , features = features
                                        , targetIndex)
                  # compute predictors for integrated model fit
                  if(integrated.model){
                    Xlist2 <- vector("list", 2)
                    for(i in 1:2){
                      Xlist2[[i]] <- Xlist[[i]]
                      tmp <- (Xlist2[[i]][1:(L-1),,drop=FALSE]+Xlist2[[i]][2:L,,drop=FALSE])/2
                      Xlist2[[i]] <- tmp*matrix(rep(diff(binned$binTime), ncol(tmp)), L-1, ncol(tmp))
                    }
                  }
                  #Data Fit B considers model on the data from all other experiments
                  parameters <- vector("list", 2)
                  for(i in 1:2){
                    dYout <- dYlist[[-i]]
                    Xout <- Xlist2[[-i]]
                    dataset <- cbind(dYout,Xout)
                    Xin <- Xlist[[i]]
                    fit <- lm(formula,data = as.data.frame(dataset))
                    #fitLM[[i]] <- fit
                    # remove coefficients resulting from singular fits (perfectly correlated predictors)
                    coefs <- coefficients(fit)
                    coefs[is.na(coefs)] <- 0
                    fitted_dY <- Xin %*% matrix(coefs, ncol(Xin), 1)
                    parameters[[i]] <- coefs

                    # Fit individual models with derivative constraint
                    fit <- constrained.smoothspline(Y[,i,drop=F]
                                                                   ,binned$binTime
                                                                   ,pen.degree=2
                                                                   ,constraint="fixed"
                                                                   ,derivative.values=fitted_dY
                                                                   ,initial.value=NA
                                                                   ,times.new=0
                                                                   ,num.folds=num.folds
                                                                   ,lambda=lambda[i])
                    RSS_B[i] <- sum(fit$residuals^2)
                    #fitSpline[[i]] <- fit
                  }
                  ## compute score
                  if(score.type=="max"){
                    if(min(abs(RSS_A))<10^-10){
                      warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type max_absolut.")
                    }
                    score <- max((RSS_B-RSS_A)/RSS_A)
                  }
                  else if(score.type=="mean"){
                    if(min(abs(RSS_A))<10^-10){
                      warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type mean_absolut.")
                    }
                    score <- mean((RSS_B-RSS_A)/RSS_A)
                  }
                  else if(score.type=="mean_absolute"){
                    score <- mean(RSS_B)
                  }
                  else if(score.type=="max_absolute"){
                    score <- max(RSS_B)
                  }
                  else if(score.type=="mean.weighted"){
                    if(min(abs(RSS_A))<10^-10){
                      warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.types mean_absolut or max_absolute.")
                    }
                    score <- mean(weight.vec*(RSS_B-RSS_A)/RSS_A)
                  }
                  else if(score.type=="stability"){
                    pairs <- combn(1:length(RSS_B), 2)
                    score <- max(apply(pairs, 2, function(x) abs(RSS_B[x[1]]-RSS_B[x[2]])))
                  }
                  else{
                    stop("Specified score.type does not exist. Use max, mean or max-mean.")
                  }
                  return(list(score = score, formula = formula))

                })
    stopCluster(cl)


    formulas <- sapply(modelScores
                       , function(singleModel){
                         return(singleModel$formula)
                       })

    modelScores <- sapply(modelScores
                          , function(singleModel){
                            return(singleModel$score)
                          })
    # scores <- vector("numeric",)
    # for (i in 1:length(modelScores)) {
    #
    # }
  }
  else{
    modelScores <- vector("numeric", num.models)
    RSS_B <- vector("numeric", 2)
    UpDown_B <- vector("numeric", 2)
    RSS3_B <- vector("numeric", 2)
    Xlist <- vector("list", 2)
    for (modelIndex in 1:num.models){
      ## collect predictors X
      index <- unique(unlist(models[[modelIndex]]))
      Xlist[[1]] <- binned$Env1[,features[index],drop=F]
      Xlist[[2]] <- binned$Env2[,features[index],drop=F]

      formula <- getFormula(model = models[[modelIndex]]
                            , features = features, targetIndex = 1)
      # compute predictors for integrated model fit
      if(integrated.model){
        Xlist2 <- vector("list", 2)
        for(i in 1:2){
          Xlist2[[i]] <- Xlist[[i]]
          tmp <- (Xlist2[[i]][1:(L-1),,drop=FALSE]+Xlist2[[i]][2:L,,drop=FALSE])/2
          Xlist2[[i]] <- tmp*matrix(rep(diff(binned$binTime), ncol(tmp)), L-1, ncol(tmp))
        }
      }
      #Data Fit B considers model on the data from all other experiments
      parameters <- vector("list", 2)
      for(i in 1:2){
        dYout <- dYlist[[-i]]
        Xout <- Xlist2[[-i]]
        dataset <- cbind(dYout,Xout)
        Xin <- Xlist[[i]]

        fit <- lm(formula,data = as.data.frame(dataset))
        # remove coefficients resulting from singular fits (perfectly correlated predictors)
        coefs <- coefficients(fit)
        coefs[is.na(coefs)] <- 0
        fitted_dY <- Xin %*% matrix(coefs, ncol(Xin), 1)
        parameters[[i]] <- coefs

        # Fit individual models with derivative constraint
        fit <- constrained.smoothspline(Y[,i,drop=F]
                                                       ,binned$binTime
                                                       ,pen.degree=2
                                                       ,constraint="fixed"
                                                       ,derivative.values=fitted_dY
                                                       ,initial.value=NA
                                                       ,times.new=0
                                                       ,num.folds=num.folds
                                                       ,lambda=lambda[i])
        RSS_B[i] <- sum(fit$residuals^2)
        UpDown_B[i] <- fit$smooth.vals[length(fit$smooth.vals)]
        RSS3_B[i] <- sum(fit$residuals[c(1, floor(L/2), L)]^2)
      }
      ## compute score
      if(score.type=="max"){
        if(min(abs(RSS_A))<10^-10){
          warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type max_absolut.")
        }
        modelScores[modelIndex] <- max((RSS_B-RSS_A)/RSS_A)
      }
      else if(score.type=="mean"){
        if(min(abs(RSS_A))<10^-10){
          warning("RSS of unconstrained smoother is very small (<10^-10). Using a relative score will lead to wrong results. Consider using the score.type mean_absolut.")
        }
        modelScores[modelIndex] <- mean((RSS_B-RSS_A)/RSS_A)
      }
      else if(score.type=="mean_absolute"){
        modelScores[modelIndex] <- mean(RSS_B)
      }
      else if(score.type=="max_absolute"){
        modelScores[modelIndex] <- max(RSS_B)
      }
      else{
        stop("Specified score.type does not exist. Use max, mean or max-mean.")
      }
    }
  }
  # modelFits <- lapply(modelScores
  #                     , function(singleModel){
  #                       return(list(lmFit = singleModel$fitLM
  #                                   , splineFit = singleModel$fitSpline
  #                                   ))
  #                     }
  #                     )
  return(list(modelScores = modelScores, formulas = formulas ))
  #return(modelScores = modelScores)
}


getFormula <- function(model, features, targetIndex){
  f <- paste(paste0("d_",features[targetIndex]),"~ ")
  for (i in 1:length(model)) {
    if(i>1)
      f <- paste(f,"+")
    f <- paste(f,paste(features[model[[i]]],collapse = ":"))
  }
  f <- paste(f,"- 1")
  return(f)
}

