DCDKM.GetModels <- function(nPredictors, option = 3){
  vars <- 1:nPredictors

  models.1 <- apply(combn(vars,2), MARGIN = 2
                    , function(x){
                      list(x[1],x[2])
                    })

  models.2 <- apply(combn(vars,2), MARGIN = 2
                    , function(x){
                      list(x)
                    })

  models.3 <- apply(combn(vars,2), MARGIN = 2
                    , function(x){
                      list(x[1],x[2],x)
                    })

  return(list(models.1 = models.1
              , models.2 = models.2
              , models.3 = models.3
  ))

}
