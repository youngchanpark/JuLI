
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

selectfun=function(cpo,bedin){
  if(grepl(';',cpo)){
    cpos=unlist(strsplit(cpo,';'))
    rnofun=function(cpos){
      cposs=unlist(strsplit(cpos,'/'))
      bedin$splno[bedin$chr==cposs[1] & bedin$pos==cposs[2] & bedin$ori==cposs[3]]
    }
    rno=unlist(sapply(cpos,rnofun))
    out=cpos[rno==max(rno)][1]
  }else{out=cpo}
  return(out)
}

