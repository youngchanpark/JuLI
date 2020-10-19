#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

NucDivFun=function(seq,conseq){
  seq=seq[order(nchar(seq))] %>% .[.!=""]
  SeqFreq=table(seq)/length(seq)
  sapply(unique(seq),function(x){
    sum(unlist(strsplit(x,"")) != unlist(strsplit(conseq,""))[1:nchar(x)])/nchar(x)*SeqFreq[names(SeqFreq)==x]
  }) %>% sum()*length(seq)/ifelse(length(seq)==1,1,(length(seq)-1))
}
