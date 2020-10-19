#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

BEDMergeFun=function(DT){
  DT=DT[,group := {ir<-IRanges(str, end); subjectHits(findOverlaps(ir, reduce(ir)))}]
  
  DT=lapply(unique(DT$group),function(x){
    DTtmp=DT[DT$group==x]
    data.frame(chr=DTtmp$chr[1],str=min(DTtmp$str),end=max(DTtmp$end),group=x,no=sum(DT$group==x),stringsAsFactors=F)}) %>% rbindlist()
  return(DT)
}