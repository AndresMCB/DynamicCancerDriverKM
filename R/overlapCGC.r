overlapCGC<-function(geneScore = NULL
                          , top = c(50,100,150,200,250))
  {
  aux <- list()
  features <-geneScore$features
  aux$CGC.baseline = intersect(geneScore$features ,CGC.driverNames$Ensembl.ID)
  geneScore <- geneScore %>%
    arrange(desc(score))%>%
    filter(score>0)


  summaryTable <- data.frame(CGC.baseline = length(aux$CGC.baseline))

  for (i in top) {
    varName <- paste0("CGC_in_top_", i)
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
