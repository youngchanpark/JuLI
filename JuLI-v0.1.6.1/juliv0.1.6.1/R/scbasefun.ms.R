#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

scbasefun.ms=function(cigar,base){
  slen=suppressWarnings(as.numeric(sub('.+M([0-9]+)S$', '\\1', cigar)))
  sbase=substring(base, nchar(base) - slen + 1)
  sbase[is.na(sbase)] = ""
  return(sbase)
}

