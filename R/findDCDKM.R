#' function to find Dynamic Cancer Drivers Kinetic Models.
#'
#' put some description here....
#'
#'
#'

findDCDKM <- function(GE1,GE2 = NULL, pTime1, pTime2
                      , GT=NULL, PPI = NULL
                      , maineffect.models = FALSE, perBin=20){

  if(!require(AMCBGeneUtils))
    devtools::install_github("AndresMCB/AMCBGeneUtils")

  # if not provided, we use Cosmic Cancer Gene Census (CGC) as conservative GT
  if(is.null(GT)){
    GT <- AMCBGeneUtils::CGCv96
  }


}
