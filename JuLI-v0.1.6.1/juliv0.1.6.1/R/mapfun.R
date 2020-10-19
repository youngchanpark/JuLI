#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

mapfun=function(gn,info,bp,gpmap,pfam,ref){
  pref=ref[ref$V2==gsub('_[[:alnum:]]+[(p].+$','',info) & ref$V13==gn,][1,]
  oes=as.numeric(unlist(strsplit(pref$V10,",")))
  oee=as.numeric(unlist(strsplit(pref$V11,",")))
  st=en=Domain=rep(NA,length(oes))
  if(pref$V4=='+'){
    exon=oee-oes
    nexon=exon/sum(exon)*0.2
    intron=oes[-1]-oee[-length(oee)]
    nintron=intron/sum(intron)*0.2
    for(n in 1:length(exon)){
      if(n==1){
        st[n]=0.05
        en[n]=st[n]+nexon[n]
      }else{
        st[n]=en[n-1]+nintron[n-1]
        en[n]=st[n]+nexon[n]
      }
    }
    if(grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=en[no]+nintron[no]*((bp-(oee)[no])/intron[no])
    }
    if(!grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=st[no]+nexon[no]*((bp-(oes)[no])/exon[no])
    }
  }
  if(pref$V4=='-'){
    exon=rev(oee-oes)
    nexon=exon/sum(exon)*0.2
    intron=rev(oes[-1]-oee[-length(oee)])
    nintron=intron/sum(intron)*0.2
    for(n in 1:length(exon)){
      if(n==1){
        st[n]=0.05
        en[n]=st[n]+nexon[n]
      }else{
        st[n]=en[n-1]+nintron[n-1]
        en[n]=st[n]+nexon[n]
      }
    }
    if(grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=en[no]+nintron[no]*((rev(oee)[no]-bp)/intron[no])
    }
    if(!grepl('Intron',info)){
      no=as.numeric(gsub('^.+[(]','',gsub('/.+$','',info)))
      nbp=st[no]+nexon[no]*((rev(oes)[no]-bp)/exon[no])
    }
  }
  if(sum(gpmap[,2]==gn)!=0){
    ppfam=pfam[grepl(gpmap[gpmap[,2]==gn,4][1],pfam$V1),]
    if( nrow(ppfam)!=0){
      index=((oes <= pref$V7 & pref$V7 < oee) | (pref$V7 < oes & oee < pref$V8) | (oes < pref$V8 & pref$V8 <=oee))
      es=oes[index]
      es[1]=pref$V7
      ee=oee[index]
      ee[sum(index)]=pref$V8
      Domain=rep(NA,length(es))
      if(pref$V4=='+'){exon=ee-es}else{exon=rev(ee-es)}
      cumexon=cumsum(exon)
      for(i in 1:nrow(ppfam)){
        dstr=unlist(strsplit(gsub('[[:space:]]+.+$','',gsub('^.+/','',ppfam$V1[i])),'-'))
        if( all(!grepl('[[:alpha:]]',dstr)) ){
          dst=as.numeric(dstr[1])*3-2
          den=as.numeric(dstr[2])*3
          Domain[(dst <= c(1,cumexon[-length(cumexon)]+1) & cumexon <= den )|(c(1,cumexon[-length(cumexon)]+1)<=dst & dst<=cumexon)|(c(1,cumexon[-length(cumexon)]+1)<=den & den<=cumexon)]=gsub('^.+[[:space:]]+','',ppfam$V2[i])
        }
      }
      utr=which(index==FALSE)
      utr5=utr[utr<(which(index==TRUE)[1])]
      utr3=utr[which(index==TRUE)[sum(index)]<utr]
      Domain=c(rep(NA,length(utr5)),Domain,rep(NA,length(utr3)))
    }
  }
  out=data.frame(st,en,Domain,nbp,stringsAsFactors=F)
  return(out)
}

