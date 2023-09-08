DCDKM.GetModels <- function(nPredictors, IncludedModels = c(1,2,3)){
  vars <- 1:nPredictors

  Collab <- NULL
  Inter <- NULL
  MainEffects <- NULL
  outcome <- NULL

  # pairwise gene collaboration linear models
  if(1%in%IncludedModels){
    outcome$Collab <- apply(combn(vars,2), MARGIN = 2
                  , function(x){
                    list(x[1],x[2])
                  })
  }

  # pairwise interaction models
  if(2%in%IncludedModels){
    outcome$Inter <- apply(combn(vars,2), MARGIN = 2
                    , function(x){
                      list(x)
                    })
  }

  # Main Effect models (i.e interactions + collaborations)
  if(3%in%IncludedModels){
    outcome$MainEffects <- apply(combn(vars,2), MARGIN = 2
                    , function(x){
                      list(x[1],x[2],x)
                    })
  }
  return(outcome)

}
