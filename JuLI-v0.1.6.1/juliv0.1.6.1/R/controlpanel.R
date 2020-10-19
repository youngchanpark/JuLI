#' common_functions
#'
#' common_functions
#' @param Defaults to NULL
#' @keywords common_functions()
#' @export
#' @examples
#' common_functions()

controlpanel=function(ControlBams=c(),
                      ID=NULL,
                      OutputPath=NULL,
                      Thread=1,
                      SplitCutoff=4){
  
  registerDoMC(Thread)
  options(scipen = 999)
  SplitCutoff = ifelse(SplitCutoff < 4,4,SplitCutoff)
  
  beds=list()
  
  for( i in 1:length(ControlBams)){
    writeLines(paste("Processing number:",i,format(Sys.time(),"%Y-%b-%d %H:%M:%S")))
    bam.info=seqinfo(BamFile(ControlBams[i]))
    sl=seqlengths(bam.info)
    chr=seqnames(bam.info) %>% .[. %in% paste0('chr',c(c(1:22),'X','Y'))]
    
    beds[[i]]=foreach(c=1:length(chr),.options.multicore=list(preschedule=FALSE)) %dopar% {
      what=c('pos','cigar')
      param=ScanBamParam(which=GRanges(chr[c],IRanges(1,sl[chr[c]])),flag=scanBamFlag(isUnmappedQuery=FALSE),what=what)
      tmp=scanBam(ControlBams[i],param=param) %>% rbindlist() %>% .[grepl('[SH]',cigar),]
      
      if(nrow(tmp)!=0){
        dat=tmp %>% mutate(end=pos+sapply(cigar,endfun)-1)
        st.br=table(dat$pos[grepl('[SH][[:digit:]]+M',dat$cigar)])
        if(length(st.br)!=0){st.br=data.table(st.br) %>% mutate(ori=0)}else{st.br=data.frame(matrix(ncol=3))[0,]}
        en.br=table(dat$end[grepl('M[[:digit:]]+[SH]',dat$cigar)])
        if(length(en.br)!=0){en.br=data.table(en.br) %>% mutate(ori=1)}else{en.br=data.frame(matrix(ncol=3))[0,]}
        sh.br=rbind(st.br,en.br) %>% filter(N >= SplitCutoff) %>% mutate(chr=chr[c]) %>% select(chr,V1,ori,N)
      }else{sh.br=data.table(matrix(ncol=4))[0,]}
      return(sh.br)
    } %>% rbindlist() %>% setNames(c('chr','pos','ori','splno'))
  } 
  
  beds=beds%>% rbindlist() %>% select(chr,pos,ori) %>% unique()
  
  fwrite(beds,paste0(OutputPath,"/",ID,".controlpanel.txt"),sep='\t',showProgress=F)
}
  
  
  
  
  