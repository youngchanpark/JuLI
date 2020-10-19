
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

dupfun=function(no,bed_tmp,bedin){
  cpo=unlist(strsplit(bed_tmp$part[no],'/'))
  data.frame(bedin[bedin$chr==cpo[1] & bedin$pos==cpo[2] & bedin$ori==cpo[3],1:5],part=paste(bed_tmp$chr[no],bed_tmp$pos[no],bed_tmp$ori[no],sep='/'),stringsAsFactors = F)
}

