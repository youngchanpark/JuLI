#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

OutputGenFun=function(Output,OutputPath,TestID){
  if(nrow(Output)!=0){
    skip=c()
    event=list()
    
    for( i in 1:nrow(Output)){
      if( sum(skip==i)==0 ){
        a=data.frame(Output$chr[i],Output$pos[i],Output$ori[i],Output$dis[i],Output$spl[i],stringsAsFactors = F)
        b=unlist(strsplit(as.character(Output$counter[i]),'/'))
        b.no=which(Output$chr==b[1]&Output$pos==b[2]&Output$ori==b[3]&Output$counter==paste(Output$chr[i],Output$pos[i],Output$ori[i],sep='/'))
        if( length(b.no)!=0 & !(a[1]==b[1] & a[2]==b[2])){
          b=Output[b.no,]
          skip=c(skip,b.no)
          b=data.frame(b$chr,b$pos,b$ori,b$dis,b$spl,stringsAsFactors = F)
          c=data.frame(a,b,stringsAsFactors = F)
          colnames(c)=paste0('V',c(1:ncol(c)))
          
          if(c$V1!=c$V6){Event='Interchromosomal_translocation'}
          if(c$V1==c$V6){
            if(c$V3==c$V8){Event='Inversion'}
            if(c$V3!=c$V8){
              if( c(c$V3,c$V8)[c(c$V2,c$V7)==min(c(c$V2,c$V7))]==1 ){Event='Deletion'}else{Event='Tandem'}
            }
          }
          c=data.frame(c,Event,stringsAsFactors = F)
          
          if( length(unique(ref$V4[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]))== 2 ){
            pref=ref[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6),]
            prefp=names(sort(table(pref$V13[pref$V4=='+']),decreasing=T))[1]
            prefn=names(sort(table(pref$V13[pref$V4=='-']),decreasing=T))[1]
            GeneA=c(prefp,prefn)
            StrandA=c('+','-')
          }
          if( length(unique(ref$V4[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]))== 1 ){
            GeneA=names(sort(table(ref$V13[ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6)]),decreasing=T))[1]
            GeneA=c(GeneA,paste0('Compl_',GeneA))
            str=unique(ref$V4[GeneA[1]==ref$V13])
            StrandA=c(str,c('-','+')[str!=c('-','+')]) 
          }
          if( sum(ref$V3==c$V1 & (ref$V5 <= c$V2 & c$V2 <= ref$V6) )== 0 ){
            pref=ref[ref$V3==c$V1,]
            mn=min(abs(pref$V5 - c$V2),abs(c$V2 - pref$V6))
            bp1=pref[(abs(pref$V5 - c$V2)==mn | abs(c$V2 - pref$V6)==mn),]
            GeneA=c(paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'),paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'))
            StrandA=c('+','-')
          }
          if( length(unique(ref$V4[ref$V3==c$V6 & (ref$V5 <= c$V7 & c$V7 <= ref$V6)]))==2 ){
            pref=ref[ref$V3==c$V6 & (ref$V5 <= c$V7 & c$V7 <= ref$V6),]
            prefp=names(sort(table(pref$V13[pref$V4=='+']),decreasing=T))[1]
            prefn=names(sort(table(pref$V13[pref$V4=='-']),decreasing=T))[1]
            GeneB=c(prefp,prefn)
            StrandB=c('+','-')
          }
          if( length(unique(ref$V4[ref$V3==c$V6 & (ref$V5 <= c$V7 & c$V7 <= ref$V6)]))==1  ){
            GeneB=names(sort(table(ref$V13[ref$V3==c$V6 & (ref$V5 <= c$V7 & c$V7 <= ref$V6)]),decreasing=T))[1]
            GeneB=c(GeneB,paste0('Compl_',GeneB))
            str=unique(ref$V4[GeneB[1]==ref$V13])
            StrandB=c(str,c('-','+')[str!=c('-','+')])      
          }
          if( sum(ref$V3==c$V6 & (ref$V5 <= c$V7 & c$V7 <= ref$V6) )==0  ){
            pref=ref[ref$V3==c$V6,]
            mn=min(abs(pref$V5 - c$V7 ),abs(c$V7  - pref$V6))
            bp1=pref[(abs(pref$V5 - c$V7)==mn | abs(c$V7  - pref$V6)==mn),]
            GeneB=c(paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'),paste('Flanking',paste(unique(bp1$V13),collapse='&'),sep='_'))
            StrandB=c('+','-')
          }
          
          cg1=data.frame(c,GeneA[1],StrandA[1],GeneB[1],StrandB[1],stringsAsFactors=F)
          cg2=data.frame(c,GeneA[2],StrandA[2],GeneB[1],StrandB[1],stringsAsFactors=F)
          cg3=data.frame(c,GeneA[1],StrandA[1],GeneB[2],StrandB[2],stringsAsFactors=F)
          cg4=data.frame(c,GeneA[2],StrandA[2],GeneB[2],StrandB[2],stringsAsFactors=F)
          colnames(cg1)=colnames(cg2)=colnames(cg3)=colnames(cg4)
          event[[i]]=rbind(cg1,cg2,cg3,cg4) %>% setNames(paste0('V',c(1:ncol(.))))
        }
      }
    }
    
    dat=rbindlist(event)
    filter=rep(T,nrow(dat))
    
    for( i in 1:nrow(dat)){
      if( dat$V3[i]==dat$V8[i] & dat$V13[i]!=dat$V15[i] ){filter[i]=F}
      if( dat$V3[i]!=dat$V8[i] & dat$V13[i]==dat$V15[i] ){filter[i]=F}
    }
    dat=dat[filter==F,]
    
    filter=rep(T,nrow(dat))
    readsum=dat$V4 + dat$V5 + dat$V9 + dat$V10
    
    for( i in 1:nrow(dat)){
      index=dat$V1[i]==dat$V1 & dat$V6[i]==dat$V6 & dat$V3[i]==dat$V3 & dat$V8[i]==dat$V8 & (dat$V2[i]-(readlen/2)) < dat$V2 & dat$V2 < (dat$V2[i]+(readlen/2)) & (dat$V7[i]-(readlen/2)) < dat$V7 & dat$V7 < (dat$V7[i]+(readlen/2))
      if( any(index) ){
        if(readsum[i] != max(readsum[index])){
          filter[i]=F
        }
      }
    }
    dat=dat[filter==T,]
  }else{
    dat=data.frame(matrix(ncol=15))[0,]
  }
  
  dat=dat %>% setNames(c('ChrA','BreakA','OriA','DisA','SplitA','ChrB','BreakB','OriB','DisB','SplitB','Event','GeneA','StrGeneA','GeneB','StrGeneB'))
  fwrite(dat,paste0(OutputPath,'/',TestID,'.txt'),sep="\t",showProgress=F)
}
