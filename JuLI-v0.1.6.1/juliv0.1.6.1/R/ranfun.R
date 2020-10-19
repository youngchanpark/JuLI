
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

ranfun=function(chr,ou){
  tab.chr=table(chr)
  bed.chr=unique(chr)
  out=foreach(i=1:length(bed.chr),.options.multicore=list(preschedule=FALSE)) %dopar% {
    paste(bed.chr[i],seq(0,tab.chr[bed.chr[i]]-0.5,ou)+1,c(seq(0,tab.chr[bed.chr[i]]-0.5,ou)[-1],tab.chr[bed.chr[i]]),sep='/')
  }
  unlist(out)
}

