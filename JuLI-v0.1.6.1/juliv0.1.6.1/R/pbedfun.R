#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

pbedfun=function(bed,bed.ran,len,medis){
  pbed.ran=unlist(strsplit(bed.ran,'/'))
  pbed=bed %>% filter(chr==pbed.ran[1]) %>% arrange(pos)
  pbed=pbed[c(as.numeric(pbed.ran[2]):as.numeric(pbed.ran[3])),]
  str=end=pbed$pos
  str[pbed$ori==0]=str[pbed$ori==0]-len
  end[pbed$ori==0]=end[pbed$ori==0]+2*medis
  str[pbed$ori==1]=str[pbed$ori==1]-2*medis
  end[pbed$ori==1]=end[pbed$ori==1]+len
  str[str<0]=end[end<0]=1
  out=data.table(chr=pbed.ran[1],str,end,pbed[,-1])
  return(out)
}
