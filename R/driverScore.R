#' Function for scoring genes in terms of their appearance in invariant models


driverScore <- function(models, features){
  d <- length(features)
  gScore <- data.frame(features = features, score = numeric(d))
  for (i in 1:d) {
    Vi <- sapply(models
                 , function(model,gene){
                   return(any(gene %in% unlist(model)))
                 }
                 ,i)
    gScore[i,2] <- sum(Vi)/d
  }
  return(gScore)
}
