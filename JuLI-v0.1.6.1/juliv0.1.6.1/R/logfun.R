#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

logfun=function(LogVal,bedLogform,bedform,col,txt){
  if(LogVal){
    idx = !(paste(bedLogform$chr,bedLogform$pos,bedLogform$ori) %in% paste(bedform$chr,bedform$pos,bedform$ori)) & bedLogform[,col]=="Pass"
    bedLogform[idx,col] = txt
  }
  bedLogform
}
