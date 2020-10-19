#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

IndexFun=function(CaseBam,ControlBam){
  
  if(BamFile(CaseBam)$index %>% isEmpty()){
    indexBam(CaseBam) %>% invisible()
  }
  if(!isEmpty(ControlBam)){
    if(BamFile(ControlBam)$index %>% isEmpty()){
      indexBam(ControlBam) %>% invisible()
    }
  }
}

