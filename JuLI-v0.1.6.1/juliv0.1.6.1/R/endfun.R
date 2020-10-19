
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

endfun=function(cigar){
  sum(as.numeric(unlist(strsplit(gsub('[[:digit:]]+[SHIN]','',cigar),'[MD]'))))
}

