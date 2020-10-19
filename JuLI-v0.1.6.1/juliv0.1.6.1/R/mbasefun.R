
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

mbasefun=function(cigar,base){
  ss=suppressWarnings(as.numeric(sub('S[[:alnum:]]+$','',cigar)))
  ss[is.na(ss)] = 0
  se=suppressWarnings(as.numeric(sub('.+M([0-9]+)S$', '\\1', cigar)))
  se[is.na(se)] = 0
  substr(base, ss + 1, nchar(base) - se)
}

