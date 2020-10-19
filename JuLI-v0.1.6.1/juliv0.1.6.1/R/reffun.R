
#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

reffun=function(gn,pos,ref){
  if( !grepl('Flanking',gn) & !grepl('Compl',gn) ){
    index= ref$V13==gn & ref$V5 <= pos & pos <= ref$V6
    gn_len=as.numeric(ref$V6[index])-as.numeric(ref$V5[index])
    if( any(grepl('NM',ref$V2[index])) ){
      gnnm=ref$V2[index][gn_len==max(gn_len[grepl('NM',ref$V2[index])]) & grepl('NM',ref$V2[index])]
      gnnm=gnnm[grepl((gsub('^.+_','',gnnm) %>% as.numeric() %>% .[.==min(.)]),gnnm)][1]
    }else{
      gnnm=ref$V2[index][gn_len==max(gn_len)]
      gnnm=gnnm[grepl((gsub('^.+_','',gnnm) %>% as.numeric() %>% .[.==min(.)]),gnnm)][1]
    }
    gnref=ref[ref$V2==gnnm & ref$V5 <= pos & pos <= ref$V6, ][1,]    
    es=as.numeric(unlist(strsplit(gnref$V10,',')))
    ee=as.numeric(unlist(strsplit(gnref$V11,',')))
    is=ee[-length(ee)]-1
    ie=es[-1]+1
    cod=as.numeric(unlist(strsplit(gnref$V16,',')))
    if( sum(es <= pos & pos <= ee)==1 ){
      if( gnref$V4 == '+' ){
        gnfun=paste0('_Exon(',which(es <= pos & pos <= ee),'/',gnref$V9,')')
        cod=paste0('_Frame(',cod[which(es <= pos & pos <= ee)],',',
                   cod[which(es <= pos & pos <= ee)+1],')')
      }
      if( gnref$V4 == '-' ){
        gnfun=paste0('_Exon(',gnref$V9-which(es <= pos & pos <= ee)+1,'/',gnref$V9,')')
        cod=paste0('_Frame(',rev(cod)[gnref$V9-which(es <= pos & pos <= ee)+1],',',
                   rev(cod)[gnref$V9-which(es <= pos & pos <= ee)],')')
      }
    }
    if( sum(is <= pos & pos <= ie)==1 ){
      if( gnref$V4 == '+' ){
        gnfun=paste0('_Intron(',which(is <= pos & pos <= ie),'/',gnref$V9-1,')')
        cod=paste0('_Frame(',cod[which(is <= pos & pos <= ie)+1],')')
      }
      if( gnref$V4 == '-' ){
        gnfun=paste0('_Intron(',gnref$V9-which(is <= pos & pos <= ie),'/',gnref$V9-1,')')
        cod=paste0('_Frame(',rev(cod)[gnref$V9-which(is <= pos & pos <= ie)+1],')')
      }
    }
    if( gnref$V5 <= pos & pos < gnref$V7 & grepl('Exon',gnfun)){
      if( gnref$V4 == '+' ){gnfun=gsub('_Exon','_5pUTR',gnfun)}
      if( gnref$V4 == '-' ){gnfun=gsub('_Exon','_3pUTR',gnfun)}
      cod=''
    }
    if( gnref$V8 < pos & pos <= gnref$V6 & grepl('Exon',gnfun)){
      if( gnref$V4 == '+' ){gnfun=gsub('_Exon','_3pUTR',gnfun)}
      if( gnref$V4 == '-' ){gnfun=gsub('_Exon','_5pUTR',gnfun)}
      cod=''
    }
    if( gnref$V8 < pos & pos <= gnref$V6 & grepl('Exon',gnfun)){
      if( gnref$V4 == '+' ){gnfun=gsub('_Exon','_3pUTR',gnfun)}
      if( gnref$V4 == '-' ){gnfun=gsub('_Exon','_5pUTR',gnfun)}
      cod=''
    }
    if(!grepl('NM',gnnm)){
      gnfun=cod=""
    }
    info=paste0(gnnm,gnfun,cod)
  }else{
    info=NA
  }
  return(info)
}
