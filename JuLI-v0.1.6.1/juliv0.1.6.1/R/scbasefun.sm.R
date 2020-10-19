#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

scbasefun.sm=function(cigar,base){
  slen=suppressWarnings(as.numeric(sub('S[[:alnum:]]+$','',cigar)))
  sbase=substr(base, 1, slen) %>% reverse()
  sbase[sbase=="NA"] = ""
  return(sbase)
}

