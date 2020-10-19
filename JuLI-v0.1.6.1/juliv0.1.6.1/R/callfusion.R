#' A function for fusion detection v0.1.6
#'
#' Fusion detection.
#' @param Defaults to NULL
#' @keywords callfusion()
#' @export
#' @examples
#' callfusion()

callfusion=function(CaseBam=NULL,
                    ControlBam=NULL,
                    TestID=NULL,
                    OutputPath=NULL,
                    Thread=1,
                    ControlPanel=NULL,
                    TargetBed=NULL,
                    Refgene=NULL,
                    Gap=NULL,
                    Reference=NULL,
                    AnalysisType='DP',
                    AnalysisUnit=1000,
                    MinMappingQuality=0,
                    SplitCutoff=2,
                    DiscordantCutoff=3,
                    NucDiv=TRUE,
                    SplitRatio=0.7,
                    MatchBase=10,
                    Log=FALSE){

  if(isEmpty(CaseBam)){MessageFun(1)}
  if(any(grepl(',',CaseBam))){CaseBam=unlist(strsplit(CaseBam,','))}
  if(any(grepl(',',TestID))){TestID=unlist(strsplit(TestID,','))}
  if(isEmpty(TestID)|length(CaseBam)!=length(TestID)){TestID = gsub('^.+/','',CaseBam)}
  if(isEmpty(OutputPath)){OutputPath = gsub(paste0('/',gsub('^.+/','',CaseBam[1])),'',CaseBam[1])}

  Thread=as.numeric(Thread);AnalysisUnit=as.numeric(AnalysisUnit);MinMappingQuality=as.numeric(MinMappingQuality);SplitCutoff=as.numeric(SplitCutoff);DiscordantCutoff=as.numeric(DiscordantCutoff);SplitRatio=as.numeric(SplitRatio);MatchBase=as.numeric(MatchBase)
  NucDiv = as.logical(NucDiv);Log = as.logical(Log)
  registerDoMC(Thread)
  options(scipen = 999)
  SplitCutoff = ifelse(SplitCutoff < 2,2,SplitCutoff)
  DiscordantCutoff = ifelse(SplitCutoff < 3,3,DiscordantCutoff)
  mat=matrix(0,5,5,dimnames=list(c("A","T","G","C","N"),c("A","T","G","C","N")))
  mat["A","A"]=mat["C","C"]=mat["G","G"]=mat["T","T"]=1

  if(isEmpty(Refgene)|isEmpty(Gap)){MessageFun(2)}
  ref <<- fread(Refgene,showProgress = F) %>% setNames(paste0('V',c(1:ncol(.))))
  GapData <<- fread(Gap,showProgress = F) %>% .[,c(2,3,4)] %>% setNames(c('chr','str','end'))
  if(isEmpty(Reference)){MessageFun(6)}
  in.fa <<- FaFile(Reference)

  invisible(lapply(CaseBam,function(x){IndexFun(x,ControlBam)}))

  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ",
                    paste0('callfusion(','CaseBam=',paste(CaseBam,collapse=','),',ControlBam=',ifelse(isEmpty(ControlBam),'NULL',ControlBam),',TestID=',paste(TestID,collapse=','),',OutputPath=',OutputPath,
                           ',Thread=',Thread,',ControlPanel=',ifelse(isEmpty(ControlPanel),'NULL',ControlPanel),',TargetBed=',ifelse(isEmpty(TargetBed),'NULL',TargetBed),
                           ',Refgene=',Refgene,',Gap=',Gap,',Reference=',Reference,',AnalysisType=',AnalysisType,',AnalysisUnit=',AnalysisUnit,',MinMappingQuality=',MinMappingQuality,
                           ',SplitCutoff=',SplitCutoff,',DiscordantCutoff=',DiscordantCutoff,',NucDiv=',as.character(NucDiv),',SplitRatio=',SplitRatio,',MatchBase=',MatchBase,',Log=',as.character(Log),')')))

  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","JuLI v0.1.6"))

  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Measuring bam statistics"))

  Stats <<- list()
  for(m in seq_len(length(CaseBam))){
    InputBam=CaseBam[m]
    bam.info=seqinfo(BamFile(InputBam))
    sl=seqlengths(bam.info)
    chr=seqnames(bam.info) %>% .[. %in% c(paste0('chr',c(c(1:22),'X','Y')),c(c(1:22),'X','Y'))] #update v0.1.4: adding chromosome name style

    readlen<<-scanBam(BamFile(InputBam,yieldSize=1000),param=ScanBamParam(what='seq')) %>% data.frame() %>% .[,1] %>% nchar() %>% median()

    out=foreach(c=1:length(chr),.options.multicore=list(preschedule=FALSE)) %dopar% {
      tmp=scanBam(InputBam,param=ScanBamParam(which=GRanges(chr[c],IRanges(1,sl[chr[c]])),what=c('cigar','isize'),mapqFilter=MinMappingQuality)) %>% rbindlist()
      isize=tmp$isize %>% abs()
      TotalReadNo=nrow(tmp)
      SplitReadNo=tmp$cigar %>% .[.!='100M'] %>% .[grepl('S',.)] %>% length()
      return(list(isize,TotalReadNo,SplitReadNo))
    }

    medis <<- out %>% sapply(.,function(x){x[[1]]}) %>% unlist() %>% median(.,na.rm=T)
    TotalReadNumber= out %>% sapply(.,function(x){x[[2]]}) %>% sum()
    SplitReadNumber= out %>% sapply(.,function(x){x[[3]]}) %>% sum()

    StatDat=data.table(Chromosome=paste(chr,collapse=';'),ReferenceLength=sum(as.numeric(sl[names(sl) %in% chr])),TotalReadNumber,SplitReadNumber,MedianInsertSize=medis,ReadLength=readlen)
    fwrite(StatDat,paste0(OutputPath,"/",TestID[m],".BamStat.txt"),sep='\t',showProgress=F)

    Stats[[m]]=list(chr=chr,readlen=readlen,medis=medis)
  }

  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Identifying candidate breaks"))

  chr = lapply(seq_len(length(CaseBam)),function(x){Stats[[x]]$chr}) %>% unlist() %>% unique()

  if(!any(grepl('chr',chr))){  #update v0.1.4
    GapData$chr=gsub('chr','',GapData$chr)
  }

  GapData %>% setkey(chr,str,end)

  if(!isEmpty(TargetBed)){
    TargetBedDat=fread(TargetBed,sep='\t',showProgress=F) %>% setNames(c('chr','str','end')) %>% setkey(chr)
    chr=chr[chr %in% unique(TargetBedDat$chr)] %>% as.character(.)
  }

  bedsplno=foreach(c = seq_len(length(chr)),.options.multicore=list(preschedule=FALSE)) %dopar% {
    what=c('pos','cigar',"mrnm","mpos")
    sh.brs=list()
    for(m in seq_len(length(CaseBam))){
      sh.br=data.table(matrix(ncol=4))[0,]
      if( any(Stats[[m]]$chr == chr[c]) ){

        InputBam = CaseBam[m]

        if(isEmpty(TargetBed)){
          param=ScanBamParam(which=GRanges(chr[c],IRanges(1,sl[chr[c]])),flag=scanBamFlag(isUnmappedQuery=FALSE),what=what,mapqFilter=MinMappingQuality)
          tmp=scanBam(InputBam,param=param) %>% rbindlist() %>% .[grepl('[SH]',cigar),] %>% unique()
        }else{
          TargetChr=chr[c]
          pbed=TargetBedDat[TargetChr] %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
          param=ScanBamParam(which=GRanges(pbed),flag=scanBamFlag(isUnmappedQuery=FALSE),what=what,mapqFilter=MinMappingQuality)
          tmp=scanBam(InputBam,param=param) %>% rbindlist() %>% .[grepl('[SH]',cigar),] %>% unique()
        }

        if(nrow(tmp)!=0){
          dat=tmp %>% mutate(end=pos+sapply(cigar,endfun)-1)
          st.br=en.br=data.frame(matrix(ncol=3))[0,]
          st.br=table(dat$pos[grepl('[SH][[:digit:]]+M',dat$cigar)]) %>% .[names(.) %in% unique(dat$pos[grepl('[S][[:digit:]]+M',dat$cigar)])]
          if(length(st.br)!=0){
            st.br=data.frame(V1=names(st.br),N=as.numeric(st.br),ori=0,stringsAsFactors = F)
          }
          en.br=table(dat$end[grepl('M[[:digit:]]+[SH]',dat$cigar)]) %>% .[names(.) %in% unique(dat$end[grepl('M[[:digit:]]+[S]',dat$cigar)])]
          if(length(en.br)!=0){
            en.br=data.frame(V1=names(en.br),N=as.numeric(en.br),ori=1,stringsAsFactors = F)
          }
          sh.br=rbind(st.br,en.br) %>% mutate(str=as.numeric(V1),chr=chr[c]) %>% arrange(str) %>% mutate(end=str) %>% select(chr,str,end,ori,N) %>% data.table() %>% foverlaps(.,GapData,type="within") %>% filter(is.na(str)) %>% select(chr,i.str,ori,N)
        }
      }
      sh.brs[[m]]=sh.br %>% setNames(c('chr','pos','ori','splno'))
    }

    tsh.brs=rbindlist(sh.brs) %>% select(chr,pos,ori) %>% unique()

    for(m in seq_len(length(sh.brs))){
      if( nrow(tsh.brs) !=0 ){
        if( nrow(sh.brs[[m]])!=0 ){
          tsh.brs=merge(tsh.brs,sh.brs[[m]],by=c("chr","pos",'ori'),all=T)
        }else{
          tsh.brs=data.table(tsh.brs,splno=0)
        }
      }else{
        tsh.brs=merge(tsh.brs,sh.brs[[m]],by=c("chr","pos",'ori'),all=T)
      }
      colnames(tsh.brs)[ncol(tsh.brs)]=paste0('splno_',m)
    }

    tsh.brs[is.na(tsh.brs)]=0

    sh.br=data.table(tsh.brs,splno=sapply(seq_len(nrow(tsh.brs)),function(x){sum(tsh.brs[x,-c(1:3)])}))
    return(sh.br)
  } %>% rbindlist() %>% filter(splno >= SplitCutoff)

  bed=bedsplno %>% select(chr,pos,ori,splno)

  bed$splno=unlist(bed$splno) #update for data.table 1.12.8

  bedLog=bed[,c(1,2,3)] %>% mutate(DiscordantPairLog="Pass",ProperPairLog="Pass")

  if(!isEmpty(ControlBam)){

    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Filtering breaks in control bam"))

    out=foreach(c=1:length(chr),.options.multicore=list(preschedule=FALSE)) %dopar% {

      if(isEmpty(TargetBed)){
        what=c('pos','cigar')
        param=ScanBamParam(which=GRanges(chr[c],IRanges(1,sl[chr[c]])),flag=scanBamFlag(isUnmappedQuery=FALSE),what=what,mapqFilter=MinMappingQuality)
        tmp=scanBam(ControlBam,param=param) %>% rbindlist() %>% .[grepl('[SH]',cigar),]
      }else{
        TargetChr=chr[c]
        pbed=TargetBedDat[TargetChr] %>% makeGRangesFromDataFrame(.,keep.extra.columns=TRUE,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
        what=c('pos','cigar')
        param=ScanBamParam(which=GRanges(pbed),flag=scanBamFlag(isUnmappedQuery=FALSE),what=what,mapqFilter=MinMappingQuality)
        tmp=scanBam(ControlBam,param=param) %>% rbindlist() %>% .[grepl('[SH]',cigar),]
      }

      if(nrow(tmp)!=0){
        dat=tmp %>% mutate(end=pos+sapply(cigar,endfun)-1)
        st.br=en.br=data.frame(matrix(ncol=3))[0,]
        st.br=table(dat$pos[grepl('[SH][[:digit:]]+M',dat$cigar)])
        if(length(st.br)!=0){
          st.br=data.frame(V1=names(st.br),N=as.numeric(st.br),ori=0,stringsAsFactors = F)
        }
        en.br=table(dat$end[grepl('M[[:digit:]]+[SH]',dat$cigar)])
        if(length(en.br)!=0){
          en.br=data.frame(V1=names(en.br),N=as.numeric(en.br),ori=1,stringsAsFactors = F)
        }
        sh.br=rbind(st.br,en.br) %>% filter(N >= SplitCutoff*2) %>% mutate(str=as.numeric(V1),chr=chr[c]) %>% arrange(str) %>% mutate(end=str) %>% select(chr,str,end,ori,N) %>% data.table() %>% foverlaps(.,GapData,type="within") %>% filter(is.na(str)) %>% select(chr,i.str,ori,N)
      }else{sh.br=data.table(matrix(ncol=4))[0,]}
      return(sh.br)
    } %>% rbindlist(use.names=FALSE) %>% setNames(c('chr','pos','ori','splno'))

    bed.n=out %>% select(chr,pos) %>% mutate(end=pos) %>% data.table() %>% setkey(chr,pos,end)
    bed=bed %>% mutate(end=pos) %>% select(chr,pos,end,ori,splno) %>% data.table() %>% foverlaps(.,bed.n,type="within") %>% filter(is.na(pos)) %>% select(chr,i.pos,ori,splno) %>% setNames(c('chr','pos','ori','splno'))

    bedLog=logfun(Log,bedLog,bed,4,"Break in the control")
    bedLog=logfun(Log,bedLog,bed,5,"Break in the control")
  }

  if(!isEmpty(ControlPanel)){

    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Filtering breaks in control panel"))

    ConPanbed=fread(ControlPanel,sep='\t',showProgress=F)
    bed.n=ConPanbed %>% select(chr,pos) %>% mutate(end=pos) %>% data.table() %>% setkey(chr,pos,end)
    bed=bed %>% mutate(end=pos) %>% select(chr,pos,end,ori,splno) %>% data.table() %>% foverlaps(.,bed.n,type="within") %>% filter(is.na(pos)) %>% select(chr,i.pos,ori,splno) %>% setNames(c('chr','pos','ori','splno'))

    bedLog=logfun(Log,bedLog,bed,4,"Break in the control panel")
    bedLog=logfun(Log,bedLog,bed,5,"Break in the control panel")
  }

  bedD=bed
  bedP=bed[bed$splno >= ifelse(SplitCutoff==2,3,SplitCutoff),]

  bedD=if(!grepl('D',AnalysisType)){
    bedD=bedD[0,]
  }else{
    bedD}
  bedP=if(!grepl('P',AnalysisType)){
    bedP=bedP[0,]
  }else{
    bedP}

  bed_D=bed_P=data.table(matrix(ncol=6))[0,]


  if(nrow(bedD)!= 0){
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Analysing discordant pairs"))
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Counting supporting reads"))
    bed.ran=ranfun(bedD$chr,AnalysisUnit)

    bedDdisno=foreach(b=1:length(bed.ran)) %dopar% {
      bchr=gsub('/.+$','',bed.ran[b])
      pbed=pbedfun(bedD,bed.ran[b],readlen,medis)
      pbedMerge= BEDMergeFun(pbed) %>% makeGRangesFromDataFrame(.,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
      pbed2=pbedfun(bedD,bed.ran[b],2,medis)

      what=c('qname',"pos","cigar","mrnm","mpos")
      param1=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isProperPair=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isMinusStrand=FALSE),what=what,mapqFilter=MinMappingQuality)
      param2=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isProperPair=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isMinusStrand=TRUE),what=what,mapqFilter=MinMappingQuality)
      param3=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isSecondaryAlignment=TRUE),what=what,mapqFilter=MinMappingQuality)
      param4=ScanBamParam(which=GRanges(pbedMerge),what=what,tag="SA",mapqFilter=MinMappingQuality)

      tdats=list()
      for(m in seq_len(length(CaseBam))){
        ppbed=data.table(matrix(ncol=6))[0,]

        InputBam = CaseBam[m]

        tmp1=scanBam(InputBam,param=param1) %>% rbindlist() %>% mutate(ori=1)
        tmp2=scanBam(InputBam,param=param2) %>% rbindlist() %>% mutate(ori=0)
        tmp3=scanBam(InputBam,param=param3) %>% rbindlist() %>% mutate(ori=2)
        dat=rbind(tmp1,tmp2,tmp3) %>% mutate_if(is.factor, as.character)

        tag.dat=scanBam(InputBam,param=param4) %>% .[sapply(., function(x){!isEmpty(x$tag$SA)})] %>% lapply(.,function(x){data.frame(x) %>% mutate_if(is.factor, as.character)}) %>% rbindlist()
        if(nrow(tag.dat)!=0){
          tag.dat= tag.dat %>% .[!is.na(.$SA),] %>% select(qname,pos,cigar,SA) %>% mutate(mrnm=gsub(',.+$','',SA),mpos=as.numeric(gsub(',.+$','',gsub('^[[:alnum:]]+,','',SA))),ori=2) %>% select(-SA)
        }

        tdat=rbind(dat,tag.dat)

        if(nrow(tdat)!=0){
          tdat$ori[grepl('M[[:digit:]]+[SH]',tdat$cigar)]=1
          tdat$ori[grepl('[SH][[:digit:]]+M',tdat$cigar)]=0
          tdat$ori[grepl('[SH].+M.+[SH]',tdat$cigar)]=2
          tdats[[m]]=data.table(CaseBam=m,unique(tdat))
        }
      }

      tdat=rbindlist(tdats)

      ppbed=data.table(matrix(ncol=6+length(CaseBam)))[0,]

      if(nrow(tdat)!=0){
        tdat=tdat %>% mutate(end=pos+sapply(cigar,endfun)-1)

        tdat.0=tdat[tdat$ori!=1,] %>% mutate(str2=pos,end2=pos) %>% data.table()
        pbed2.0=pbed2[pbed2$ori==0,] %>% select(chr,str,end,pos) %>% setNames(c("chr","str2","end2","pos")) %>% setkey(str2,end2)
        tdat.pbed.0=foverlaps(tdat.0,pbed2.0,type="within") %>% filter(!is.na(chr)) %>% data.table() %>% select(pos,CaseBam,qname,i.pos,cigar,mrnm,mpos,ori,end) %>% setNames(c('pos','CaseBam','qname','str','cigar','mrnm','mpos','ori','end')) %>% setkey(pos)

        tdat.1=tdat[tdat$ori!=0,] %>% mutate(str2=end,end2=end) %>% data.table()
        pbed2.1=pbed2[pbed2$ori==1,] %>% select(chr,str,end,pos) %>% setNames(c("chr","str2","end2","pos")) %>% setkey(str2,end2)
        tdat.pbed.1=foverlaps(tdat.1,pbed2.1,type="within") %>% filter(!is.na(chr)) %>% data.table() %>% select(pos,CaseBam,qname,i.pos,cigar,mrnm,mpos,ori,end) %>% setNames(c('pos','CaseBam','qname','str','cigar','mrnm','mpos','ori','end')) %>% setkey(pos)

        namefun=function(no){
          bpos=pbed$pos[no]
          bori=pbed$ori[no]
          if(bori==0){pdat=tdat.pbed.0[J(bpos)] %>% .[!(grepl('[SH]',cigar) & bpos < str)]}
          if(bori==1){pdat=tdat.pbed.1[J(bpos)] %>% .[!(grepl('[SH]',cigar) & end < bpos)]}
          pdat=pdat[!(mrnm==bchr & mpos %between% c(bpos-medis*2,bpos+medis*2))]
          pdat=pdat[order(match(pdat$qname,names(sort(table(pdat$qname),decreasing=T)))),] %>% .[!(.$qname %in% .$qname[duplicated(paste(CaseBam,str,mrnm,mpos,end))])]
          tc=intersect(table(pdat$mrnm) %>% .[.>=DiscordantCutoff] %>% names(),table(pdat$mrnm[grepl('[SH]',pdat$cigar)]) %>% names()) #update0828

          namefunout=data.table(matrix(ncol=2+length(CaseBam)))

          if(length(tc)!=0){
            pdat=pdat[pdat$mrnm %in% tc,]  %>% .[order(.$mrnm,.$mpos),]
            overfun=function(tc){
              out=NULL
              cpdat=pdat[pdat$mrnm==tc,]
              clu.pos=which(abs(c(cpdat$mpos[-1],-(medis*2))-cpdat$mpos) > medis*2)
              clunamefun=function(cp){
                out=NULL
                tmp.cpdat=cpdat[c((c(0,clu.pos[-length(clu.pos)])[clu.pos==cp]+1):clu.pos[clu.pos==cp]),]
                if( length(unique(tmp.cpdat$qname)) >= DiscordantCutoff & any(grepl('[SH]',tmp.cpdat$cigar)) ){ ##update0828
                  uniq.qname=unique(tmp.cpdat$qname)
                  dfout=lapply(seq_len(length(CaseBam)),function(x){length(unique(tmp.cpdat$qname[tmp.cpdat$CaseBam==x]))}) %>% unlist()
                  out=list(length(uniq.qname),uniq.qname,dfout)
                }
                return(out)
              }
              out=lapply(clu.pos,clunamefun) %>% .[!isEmpty(.)]
              if(any(!isEmpty(out))){
                out=out[[which.max(sapply(seq_len(length(out)),function(x){out[[x]][[1]]}))]]
              }
              return(out)
            }
            out=lapply(tc,overfun) %>% .[!isEmpty(.)]
            if(any(!isEmpty(out))){
              out=out[[which.max(sapply(seq_len(length(out)),function(x){out[[x]][[1]]}))]]
              namefunout=data.frame(V1=out[[1]],V2=paste(out[[2]],collapse='/'),matrix(out[[3]],nrow=1),stringsAsFactors = F) %>% setNames(paste0('V',seq_len(2+length(CaseBam))))
            }
          }
          return(namefunout)
        } %>% data.frame()

        ppbed=cbind(pbed, lapply(seq_len(nrow(pbed)),namefun) %>% rbindlist()) %>% select(-str,-end,-group)
      }
      return(ppbed %>% setNames(c('chr','pos','ori','splno','disno','qnames',paste0('disno_',seq_len(length(CaseBam))))))
    } %>% rbindlist() %>% filter(disno >= DiscordantCutoff & !is.na(disno))

    bedD=bedDdisno %>% select(chr,pos,ori,splno,disno,qnames)

    bedLog=logfun(Log,bedLog,bedD,4,"Reads count under the cutoff")
  }

  if(nrow(bedD)!= 0){
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Generating consensus contigs"))
    bed.ran=ranfun(bedD$chr,AnalysisUnit)

    bedD=foreach(b=1:length(bed.ran)) %dopar% {
      bchr=gsub('/.+$','',bed.ran[b])
      pbed=pbedfun(bedD,bed.ran[b],readlen,medis)
      pbedMerge= BEDMergeFun(pbed) %>% makeGRangesFromDataFrame(.,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
      param=ScanBamParam(which=GRanges(pbedMerge),what=c("pos","cigar","mapq",'seq'),mapqFilter=MinMappingQuality)
      dats=list()
      for(m in seq_len(length(CaseBam))){

        InputBam = CaseBam[m]

        dat=scanBam(InputBam,param=param) %>% lapply(.,function(x){data.frame(x) %>% mutate_if(is.factor, as.character)}) %>% rbindlist() %>% .[grepl('[SH]',cigar),]

        if(nrow(dat)!=0){
          dat=dat %>% mutate(end=pos+sapply(cigar,endfun)-1)
        }else{
          dat=data.table(matrix(ncol=5))[0,]
        }

        dats[[m]]=dat %>% setNames(c("str","mapq","cigar","seq","end"))
      }

      dat=rbindlist(dats)

      dat_SM=dat[grepl('[SH][[:digit:]]+M',dat$cigar),] %>% data.table() %>% setkey(str)
      dat_MS=dat[grepl('M[[:digit:]]+[SH]',dat$cigar),] %>% data.table() %>% setkey(end)

      mmfun=function(bpos,bori){
        if(bori==0){
          pdat=dat_SM[J(bpos),]
          scbase=scbasefun.sm(pdat$cigar,pdat$seq)
          mbase = mbasefun(pdat$cigar,pdat$seq)
        }
        if(bori==1){
          pdat=dat_MS[J(bpos),]
          scbase=scbasefun.ms(pdat$cigar,pdat$seq)
          mbase = mbasefun(pdat$cigar,pdat$seq) %>% reverse()
        }

        con_mbase=consensusString(DNAStringSet(mbase),ambiguityMap='N',threshold=0.5)
        NucDivmbase=NucDivFun(mbase,con_mbase)
        MeanMapq=mean(pdat$mapq)

        if(any(nchar(scbase) >= MatchBase)){
          con_scbase=consensusString(DNAStringSet(scbase),ambiguityMap='N',threshold=0.5)
          NucDivscbase=NucDivFun(scbase,con_scbase)
        }else{
          con_scbase='N'
          NucDivscbase=NA
        }
        return(c(paste(con_scbase,con_mbase,sep='/'),MeanMapq,NucDivscbase,NucDivmbase))
      }

      return(cbind(pbed %>% select(chr,pos,ori,splno,disno,qnames),mapply(mmfun,pbed$pos,pbed$ori) %>% t() %>% data.frame() %>% mutate_if(is.factor, as.character)))
    } %>% rbindlist() %>% setNames(c("chr","pos","ori","splno","disno","qnames","seq","MeanMapq","nucdiv_scbase","nucdiv_mbase")) %>% mutate(MeanMapq=as.numeric(MeanMapq),nucdiv_scbase=as.numeric(nucdiv_scbase),nucdiv_mbase=as.numeric(nucdiv_mbase))

    if(NucDiv){
      nucdiv_cutoff=ifelse(sum(bedD$nucdiv_mbase!=0)>=2,mean(bedD$nucdiv_mbase[bedD$nucdiv_mbase!=0])+2*sd(bedD$nucdiv_mbase[bedD$nucdiv_mbase!=0]),0)
      bedD=bedD %>% filter((nucdiv_scbase <= nucdiv_cutoff | is.na(nucdiv_scbase)) & (nucdiv_mbase <= nucdiv_cutoff))

      bedLog=logfun(Log,bedLog,bedD,4,"Break with high nucleotide diversity")
    }
    bedD=bedD %>% select(chr,pos,ori,splno,disno,seq,qnames,MeanMapq)
  }

  if(nrow(bedD)!= 0){
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Comparing contigs between breaks"))
    nbed=rbindlist(foreach(i=1:nrow(bedD)) %dopar% {data.table(bedD[i,c('chr','pos','ori','splno','disno','seq','MeanMapq')],qname=unlist(strsplit(bedD$qnames[i],'/')))})
    setkey(nbed,"qname")

    chrlist=chr

    count.brD=foreach(b=1:nrow(bedD)) %dopar% {
      count.br1=count.br2=bbedout=NULL
      pnbed=nbed[unique(unlist(strsplit(bedD$qnames[b],'/')))] %>% select(chr,pos,ori,splno,disno,seq,MeanMapq) %>% unique() %>% arrange(pos) %>% left_join(data.table(chr=chr[chr %in% .$chr]),.,by='chr') %>% .[-(1:last(which(.$chr==bedD$chr[b] & .$pos==bedD$pos[b]))),] %>% filter(!(chr==bedD$chr[b] & pos < (bedD$pos[b]+medis)))

      if(nrow(pnbed) !=0){
        basecomfun=function(no){
          bori=pnbed$ori[no]
          bseq=pnbed$seq[no]
          rbase=unlist(strsplit(bedD$seq[b],'/'))
          pbase=unlist(strsplit(bseq,'/'))
          if(bedD$ori[b] == bori){pbase=chartr("ATGC","TACG",pbase)}
          com1=pairwiseAlignment(rbase[1],substr(pbase[2],1,nchar(rbase[1])*1.5),type="local",substitutionMatrix=mat,gapOpening=-1,gapExtension=-1,scoreOnly=T)
          com2=pairwiseAlignment(pbase[1],substr(rbase[2],1,nchar(pbase[1])*1.5),type="local",substitutionMatrix=mat,gapOpening=-1,gapExtension=-1,scoreOnly=T)
          out=(com1 >= nchar(gsub('N','',rbase[1]))*SplitRatio & com1 >=MatchBase) | (com2 >=nchar(gsub('N','',pbase[1]))*SplitRatio & com2 >= MatchBase)
          return(out)
        }
        count.br1=pnbed[sapply(seq_len(nrow(pnbed)),basecomfun),]
      }

      pbed=bedD[b,]
      rbase=unlist(strsplit(pbed$seq,'/'))

      if(nchar(rbase[1]) >= MatchBase*2){
        if(pbed$ori==0){pbed=cbind(pbed,str=pbed$pos,end=pbed$pos+medis*2)}else{pbed=cbind(pbed,str=pbed$pos-medis*2,end=pbed$pos)}
        pbedMerge= makeGRangesFromDataFrame(pbed,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
        what=c('qname',"pos","cigar","mrnm","mpos")
        param=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isProperPair=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE),what=what,mapqFilter=ifelse(MinMappingQuality==0,1,MinMappingQuality))
        dats=list()
        for(m in seq_len(length(CaseBam))){
          InputBam <<- CaseBam[m]
          dat=scanBam(InputBam,param=param) %>% rbindlist() %>% mutate_if(is.factor, as.character)
          dats[[m]]= dat %>% .[!(.$mrnm==pbed$chr & (pbed$pos-medis*4) < .$mpos & .$mpos < (pbed$pos+medis*4)),] %>% .[!duplicated(paste(.$pos,.$cigar,.$mrnm,.$mpos)),] %>% data.table(.,chr=.$mrnm,str=.$mpos-medis*2,end=.$mpos+medis*2)
        }

        dat=rbindlist(dats)

        datbedmerge=BEDMergeFun(dat) %>% .[.$chr %in% chrlist] %>% .[.$no >= DiscordantCutoff*2]

        if( nrow(datbedmerge)!=0 ){

          bedinput= makeGRangesFromDataFrame(datbedmerge,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
          what=c('qname','rname',"pos","cigar","mapq")
          param1=ScanBamParam(which=GRanges(bedinput),flag=scanBamFlag(isProperPair=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isMinusStrand=FALSE),what=what,mapqFilter=ifelse(MinMappingQuality==0,1,MinMappingQuality))
          param2=ScanBamParam(which=GRanges(bedinput),flag=scanBamFlag(isProperPair=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE,isMinusStrand=TRUE),what=what,mapqFilter=ifelse(MinMappingQuality==0,1,MinMappingQuality))
          dat2s=list()

          for(m in seq_len(length(CaseBam))){
            InputBam <<- CaseBam[m]
            tmp1=scanBam(InputBam,param=param1) %>% rbindlist() %>% mutate(ori=1)
            tmp2=scanBam(InputBam,param=param2) %>% rbindlist() %>% mutate(ori=0)
            dat2s_tmp=rbind(tmp1,tmp2) %>% mutate_if(is.factor, as.character)
            if(nrow(dat2s_tmp)!=0){
              dat2s[[m]]=data.table(CaseBam=m,dat2s_tmp)
            }
          }

          dat2=rbindlist(dat2s) %>% .[.$qname %in% dat$qname[dat$group %in% datbedmerge$group],]

          if(length(unique(dat2$qname)) >= DiscordantCutoff*2){
            tbbed=data.table(datbedmerge,lapply(seq_len(nrow(datbedmerge)),function(x){
              idx=dat2$qname %in% dat$qname[dat$group==datbedmerge$group[x]]
              ttb=data.frame(pair.no=length(unique(dat2$qname[idx])),mapq=mean(dat2$mapq[idx]),stringsAsFactors=F)
              tttb=data.frame(matrix(sapply(seq_len(length(CaseBam)),function(y){length(unique(dat2$qname[idx & dat2$CaseBam==y]))}),nrow=1))
              colnames(tttb)=paste0('disno_',seq_len(length(CaseBam)))
              return(cbind(ttb,tttb))}) %>% rbindlist()) %>% .[.$pair.no >= DiscordantCutoff*2]

            if( nrow(tbbed) !=0 ){

              pacoms=list()
              for( t in seq_len(nrow(tbbed)) ){
                bbed=tbbed[t,]
                dat3=dat2[dat2$qname %in% dat$qname[dat$group==bbed$group],]
                tdat=dat3 %>% mutate(end=pos+sapply(cigar,endfun)-1)

                tdat$ori[grepl('M[[:digit:]]+[SH]',tdat$cigar)]=1
                tdat$ori[grepl('[SH][[:digit:]]+M',tdat$cigar)]=0
                tdat$ori[grepl('[SH].+M.+[SH]',tdat$cigar)]=2

                tabtt=table(tdat$ori)
                tori=names(tabtt)[which.max(tabtt)]

                bbedinput=makeGRangesFromDataFrame(bbed,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')

                pbase=getSeq(in.fa, bbedinput)[[1]]

                if(bedD$ori[b] == tori){pbase=chartr("ATGC","TACG",pbase)}
                if(tori== 1){pbase=rev(pbase)}

                pacoms[[t]]=pairwiseAlignment(rbase[1],pbase,type="local",substitutionMatrix=mat,gapOpening=-1,gapExtension=-1)
              }

              idx=which.max(sapply(seq_len(nrow(tbbed)),function(x){pacoms[[x]]@score}))[1]
              bbed=tbbed[idx]
              pacom=pacoms[[idx]]

              if( pacom@score >= nchar(gsub('N','',rbase[1]))*SplitRatio &  pacom@score >=MatchBase*2){
                if(tori==0){
                  brpos=bbed$str+pacom@subject@range@start-1
                  if( all(brpos <= tdat$pos[tdat$ori==0])){
                    count.br2=data.frame(chr=bbed$chr,pos=brpos,ori=tori,splno=0,disno=bbed$pair.no,seq=as.character(pacom@subject),MeanMapq=bbed$mapq,stringsAsFactors = F)
                    bbedout=bbed
                  }
                }
                if(tori==1){
                  brpos=bbed$end-pacom@subject@range@start+1
                  if( all(tdat$pos[tdat$ori==1] <= brpos)){
                    count.br2=data.frame(chr=bbed$chr,pos=brpos,ori=tori,splno=0,disno=bbed$pair.no,seq=as.character(pacom@subject),MeanMapq=bbed$mapq,stringsAsFactors = F)
                    bbedout=bbed
                  }
                }
              }
            }
          }
        }
      }

      count.br=rbind(count.br1,count.br2) %>% .[which.max((.$splno+.$disno)*.$MeanMapq),]

      if(!all(isEmpty(count.br))){
        count.br=paste(paste(count.br$chr,count.br$pos,count.br$ori,sep='/'),paste(count.br$disno,count.br$MeanMapq,sep='/'),sep=';')
      }else{
        count.br=NA
      }
      list(count.br,bbedout)
    }

    bed_tmp=bedD %>% select(chr,pos,ori,splno,disno,MeanMapq) %>% mutate(part=sapply(seq_len(nrow(bedD)),function(x){count.brD[[x]][[1]]})) %>% filter(!is.na(part))

    if(nrow(bed_tmp)!=0){
      cpo=paste(bedD$chr,bedD$pos,bedD$ori,sep='/')
      btcpo=paste(bed_tmp$chr,bed_tmp$pos,bed_tmp$ori,sep='/')
      gsubpart=gsub(';.+$','',bed_tmp$part)

      bed_D=foreach(b=seq_len(nrow(bed_tmp))) %dopar% {
        outtb=NULL
        idx=gsubpart[b]==cpo
        if(any(idx)){
          outtb=rbind(bed_tmp[b,],data.frame(bedD[idx,-c(6,7)],part=btcpo[b],stringsAsFactors = F))
        }else{
          idx=gsubpart[b]==gsubpart
          if( sum(idx)==1 ){
            str=unlist(strsplit(bed_tmp$part[b],'[/;]'))
            bed_add=data.frame(chr=str[1],pos=as.numeric(str[2]),ori=as.numeric(str[3]),splno=0,disno=as.numeric(str[4]),MeanMapq=as.numeric(str[5]),part=btcpo[b],stringsAsFactors = F)
            outtb=rbind(bed_tmp[b,],bed_add)
          }else{
            tmp=bed_tmp[idx,]
            tmpcpo=paste(tmp$chr,tmp$pos,tmp$ori,sep='/')
            tmp.part=gsub('^.+;','',tmp$part)
            partmax=which.max(tmp$disno * (tmp$MeanMapq+0.1) * as.numeric(gsub('/.+$','',tmp.part)) * (as.numeric(gsub('^.+/','',tmp.part))+0.01))
            if( which(tmpcpo==btcpo[b])==partmax){
              tmptmp=tmp[partmax,]
              str=unlist(strsplit(tmptmp$part,'[/;]'))
              bed_add=data.frame(chr=str[1],pos=as.numeric(str[2]),ori=as.numeric(str[3]),splno=0,disno=as.numeric(str[4]),MeanMapq=as.numeric(str[5]),part=btcpo[b],stringsAsFactors = F)
              outtb=rbind(tmptmp,bed_add)
            }
          }
        }
        if(all(outtb$MeanMapq < 1)){
          outtb=NULL
        }
        outtb
      } %>% rbindlist()

      bed_D$part=gsub(';.+$','',bed_D$part)

      bed_D=unique(bed_D)

      cpo=paste(bed_D$chr,bed_D$pos,bed_D$ori,sep='/')
      unicpo=table(cpo)
      unicpo = names(unicpo)[unicpo!=1]
      gsubpart=gsub(';.+$','',bed_D$part)

      if(length(unicpo)!=0){
        bed_uni=foreach(b=seq_len(length(unicpo))) %dopar% {
          tmp=bed_D[cpo %in% gsubpart[cpo==unicpo[b]],] %>% .[which.max((.$splno+.$disno)*.$MeanMapq),] %>% .[1,]
          rbind(bed_D[cpo==unicpo[b] & gsubpart==paste(tmp$chr,tmp$pos,tmp$ori,sep='/'),],tmp)
        } %>% rbindlist()
        bed_D=rbind(bed_D[!(cpo %in% c(unicpo,gsubpart[cpo %in% unicpo])),],bed_uni)
      }

      bed_D=bed_D[,-6] %>% unique()
    }

    bedLog=logfun(Log,bedLog,bed_D,4,"Unmatched contig sequences")
  }


  if(nrow(bedP)!= 0){
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Analysing proper pairs"))
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Counting supporting reads"))
    bed.ran=ranfun(bedP$chr,AnalysisUnit)

    bedPsplno=foreach(b=1:length(bed.ran)) %dopar% {
      bchr=gsub('/.+$','',bed.ran[b])
      pbed=pbedfun(bedP,bed.ran[b],readlen,medis)
      pbedMerge= BEDMergeFun(pbed) %>% makeGRangesFromDataFrame(.,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')

      param=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isProperPair=TRUE,isUnmappedQuery=FALSE),what=c("qname","pos","cigar","isize"),mapqFilter=MinMappingQuality)

      ppbeds=list()
      for(m in seq_len(length(CaseBam))){
        ppbed=data.table(matrix(ncol=5))[0,]

        InputBam = CaseBam[m]

        dat=scanBam(InputBam,param=param) %>% rbindlist() %>% .[grepl('S',cigar),] %>% setNames(c("qname","str","cigar","isize"))

        if(nrow(dat)!=0){
          dat=dat %>% mutate(ori=0,end=str+sapply(cigar,endfun)-1)
          dat$ori[grepl('M[[:digit:]]+[S]',dat$cigar)]=1
          dat$ori[grepl('[S][[:digit:]]+M',dat$cigar)]=0
          dat$ori[grepl('[S].+M.+[S]',dat$cigar)]=2
          dat=unique(dat)
          dat_0=dat[dat$ori!=1,] %>% data.table() %>% setkey(str)
          dat_1=dat[dat$ori!=0,] %>% data.table() %>% setkey(end)

          namefun=function(bpos,bori){
            if(bori==0){pdat=dat_0[J(bpos),]}
            if(bori==1){pdat=dat_1[J(bpos),]}
            if(any(!is.na(pdat$qname))){
              pdat=pdat[order(match(pdat$qname,names(sort(table(pdat$qname),decreasing=T)))),]
              pdat=pdat[!duplicated(paste(pdat$str,pdat$cigar,pdat$end)),]
              namefunout=length(unique(pdat$qname))
            }else{
              namefunout=0
            }
            return(namefunout)
          }

          ppbed=pbed %>% select(chr,pos,ori) %>% mutate(splno=mapply(namefun,.$pos,.$ori),disno=0)
        }
        ppbeds[[m]]=ppbed %>% setNames(c("chr","pos","ori","splno","disno"))
      }

      tsh.brs=pbed %>% select(chr,pos,ori) %>% data.frame(.,stringsAsFactors = F)

      for(m in seq_len(length(ppbeds))){
        tsh.brs=merge(tsh.brs,ppbeds[[m]],by=c("chr","pos",'ori'),all=T)
        colnames(tsh.brs)[ncol(tsh.brs)-1]=paste0('splno_',m)
        colnames(tsh.brs)[ncol(tsh.brs)]=paste0('disno_',m)
      }

      tsh.brs[is.na(tsh.brs)]=0

      splcol=which(grepl('splno_',colnames(tsh.brs)))

      ppbed=data.table(tsh.brs,splno=sapply(seq_len(nrow(tsh.brs)),function(x){sum(tsh.brs[x,splcol],na.rm=T)}),disno=0)

      return(ppbed)
    } %>% rbindlist() %>% filter(splno >= ifelse(SplitCutoff==2,3,SplitCutoff))

    bedP=bedPsplno %>% select(chr,pos,ori,splno,disno)

    bedLog=logfun(Log,bedLog,bedP,5,"Reads count under the cutoff")
  }

  if(nrow(bedP)!= 0){
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Generating consensus contigs"))
    bed.ran=ranfun(bedP$chr,AnalysisUnit)

    bedP=foreach(b=1:length(bed.ran)) %dopar% {
      bchr=gsub('/.+$','',bed.ran[b])
      pbed=pbedfun(bedP,bed.ran[b],readlen,medis)
      pbedMerge= BEDMergeFun(pbed) %>% makeGRangesFromDataFrame(.,ignore.strand=TRUE,seqnames.field='chr',start.field='str',end.field='end')
      param=ScanBamParam(which=GRanges(pbedMerge),flag=scanBamFlag(isUnmappedQuery=FALSE),what=c("pos","cigar","seq"),mapqFilter=MinMappingQuality)
      dats=list()
      for(m in seq_len(length(CaseBam))){

        InputBam = CaseBam[m]

        dat=scanBam(InputBam,param=param) %>% lapply(.,function(x){data.frame(x) %>% mutate_if(is.factor, as.character)}) %>% rbindlist() %>% .[grepl('[SH]',cigar),]

        if(nrow(dat)!=0){
          dat=dat %>% mutate(end=pos+sapply(cigar,endfun)-1)
        }else{
          dat=data.table(matrix(ncol=4))[0,]
        }

        dats[[m]]=dat %>% setNames(c("str","cigar","seq","end"))
      }

      dat=rbindlist(dats)

      dat_SM=dat[grepl('[SH][[:digit:]]+M',dat$cigar),] %>% data.table() %>% setkey(str)
      dat_MS=dat[grepl('M[[:digit:]]+[SH]',dat$cigar),] %>% data.table() %>% setkey(end)

      mmfun=function(bpos,bori){
        if(bori==0){
          pdat=dat_SM[J(bpos)]
          scbase=scbasefun.sm(pdat$cigar,pdat$seq)
          mbase = mbasefun(pdat$cigar,pdat$seq)
        }
        if(bori==1){
          pdat=dat_MS[J(bpos)]
          scbase=scbasefun.ms(pdat$cigar,pdat$seq)
          mbase = mbasefun(pdat$cigar,pdat$seq) %>% reverse()
        }

        con_mbase=consensusString(DNAStringSet(mbase),ambiguityMap='N',threshold=0.5)
        NucDivmbase=NucDivFun(mbase,con_mbase)

        if(any(nchar(scbase) >= MatchBase)){
          con_scbase=consensusString(DNAStringSet(scbase),ambiguityMap='N',threshold=0.5)
          NucDivscbase=NucDivFun(scbase,con_scbase)

          out=c(paste(con_scbase,con_mbase,sep='/'),NucDivscbase,NucDivmbase)
        }else{out=c(NA,NA,NA)}
        return(out)
      }

      return(cbind(pbed %>% select(chr,pos,ori,splno,disno),mapply(mmfun,pbed$pos,pbed$ori) %>% t() %>% data.frame() %>% mutate_if(is.factor, as.character)))
    } %>% rbindlist() %>% setNames(c("chr","pos","ori","splno","disno","seq","nucdiv_scbase","nucdiv_mbase")) %>% mutate(nucdiv_scbase=as.numeric(nucdiv_scbase),nucdiv_mbase=as.numeric(nucdiv_mbase)) %>% .[!is.na(.$nucdiv_mbase),]

    if(NucDiv){
      nucdiv_cutoff=ifelse(sum(bedP$nucdiv_mbase!=0)>=2,mean(bedP$nucdiv_mbase[bedP$nucdiv_mbase!=0])+2*sd(bedP$nucdiv_mbase[bedP$nucdiv_mbase!=0]),0)
      tbed=bedP %>% filter(nucdiv_scbase <= nucdiv_cutoff & nucdiv_mbase <= nucdiv_cutoff)

      bedLog=logfun(Log,bedLog,tbed,5,"Break with high nucleotide diversity")
    }else{
      tbed=bedP
    }

    tbed=tbed %>% select(chr,pos,ori,splno,disno,seq)
  }

  if(nrow(bedP)!= 0){
    if(nrow(tbed)!= 0){
      writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","-Comparing contigs between breaks"))
      bed_P=list()
      for(c in 1:length(chr)){
        bed=data.frame(tbed,stringsAsFactors = F) %>% .[.$chr==chr[c],] %>% mutate(str=pos,end=pos) %>% data.table()
        pbed=data.table(str=bed$pos-medis*2,end=bed$pos+medis*2,br=bed$pos) %>% setkey(str,end)
        ppbed=foverlaps(bed,pbed,type="within") %>% select(br,chr,pos,ori,splno,disno,seq) %>% data.table() %>% .[abs(.$br-.$pos) > readlen/2] %>% setkey(br)
        bed=bed[bed$pos %in% unique(ppbed$br)]

        if(nrow(bed)!=0){
          bed.ran=ranfun(bed$chr,100) %>% gsub(paste0(chr[c],'/'),'',.)
          bed_tmp=foreach(b=1:length(bed.ran)) %dopar% {

            pbed=bed[gsub('/.+$','',bed.ran[b]):gsub('^.+/','',bed.ran[b])]
            pppbed=ppbed[J(pbed$pos),allow.cartesian=TRUE] %>% setkey(br)

            pallcomfun=function(no){
              bed.can=pppbed[J(pbed$pos[no])] %>% .[.$pos > pbed$pos[no]]
              if(nrow(bed.can)!=0 ){
                basecomfun=function(bpos,bori,bseq){
                  out=NA
                  rbase=unlist(strsplit(pbed$seq[no],'/'))
                  pbase=unlist(strsplit(bseq,'/'))
                  if(pbed$ori[no] == bori){pbase=chartr("ATGC","TACG",pbase)}
                  com1=pairwiseAlignment(rbase[1],substr(pbase[2],1,nchar(rbase[1])*1.5),type="local",substitutionMatrix=mat,gapOpening=-1,gapExtension=-1,scoreOnly=T)
                  com2=pairwiseAlignment(pbase[1],substr(rbase[2],1,nchar(pbase[1])*1.5),type="local",substitutionMatrix=mat,gapOpening=-1,gapExtension=-1,scoreOnly=T)
                  if( com1 >= nchar(rbase[1])*SplitRatio & com2 >= nchar(pbase[1])*SplitRatio ){
                    out=paste(chr[c],bpos,bori,sep='/')
                  }
                  return(out)
                }
                count.br.each=unique(unlist(mapply(basecomfun,bed.can$pos,bed.can$ori,bed.can$seq)))
                count.br=paste(count.br.each[!is.na(count.br.each)],collapse=';')
              }else{count.br=NA}
              count.br
            }
            return(pbed %>% select(chr,pos,ori,splno,disno) %>% mutate(part=unlist(lapply(c(1:nrow(pbed)),pallcomfun))))
          } %>% rbindlist() %>% filter(!is.na(part) & part!='')

          if(nrow(bed_tmp)!=0){
            bed_tmp=data.table(bed_tmp[,c(1:5)],part=unlist(foreach(b=1:nrow(bed_tmp)) %dopar% {selectfun(bed_tmp$part[b],bed)}))
            bed_tmp=rbind(bed_tmp,rbindlist(foreach(b=1:nrow(bed_tmp)) %dopar% {dupfun(b,bed_tmp,bed)}))
          }
        }else{bed_tmp=data.frame(matrix(ncol=6))[0,]}
        bed_P[[c]]=bed_tmp %>% setNames(c("chr","pos","ori","splno","disno","part"))
      }
      bed_P=rbindlist(bed_P)
      bedLog=logfun(Log,bedLog,bed_P,5,"Unmatched contig sequences")
    }
  }

  if(Log){
    fwrite(bedLog,paste0(OutputPath,'/',paste0(TestID,collapse='_'),'.log'),sep="\t",showProgress=F)
  }

  colnames(bed_D)=colnames(bed_P)=c('chr','pos','ori','spl','dis','counter')

  Output = rbind(bed_D,bed_P) %>% filter(!duplicated(paste(chr,pos,ori,counter)))

  if(any(Output$spl==0)){
    addtb=Output[Output$spl==0,]
    bedDcpo=paste(bedD$chr,bedD$pos,bedD$ori,sep='/')
    addtbtb=lapply(seq_len(nrow(addtb)),function(x){count.brD[[which(bedDcpo==addtb$counter[x])]][[2]]}) %>% rbindlist()
    addtb=data.frame(addtb,addtbtb,stringsAsFactors = F)
  }else{
    addtb=NULL
  }

  if(grepl('D',AnalysisType) & nrow(bed_D)!=0){
    bedDmerge=merge(Output[Output$dis!=0 & Output$spl!=0,],bedsplno,by=c('chr','pos','ori')) %>% merge(.,bedDdisno,by=c('chr','pos','ori')) %>% data.frame(.,stringsAsFactors = F)
    bedDmerge=bedDmerge[,-which(grepl('qname',colnames(bedDmerge)))]
  }
  if(grepl('P',AnalysisType) & nrow(bed_P)!=0){
    bedPmerge=merge(Output[Output$dis==0,],bedPsplno,by=c('chr','pos','ori')) %>% data.frame(.,stringsAsFactors = F)
  }

  OutputList=list()
  for(m in seq_len(length(CaseBam))){
    bed_Deach=bed_Peach=bed_Deachadd=data.frame(matrix(ncol=6))[0,]
    if(grepl('D',AnalysisType) & nrow(bed_D)!=0){
      bed_DeachOri=bedDmerge[,c(1,2,3,which(grepl(paste0('_',m),colnames(bedDmerge))),6)]
      bed_Deach=bed_DeachOri %>% .[.[,5] >= 1,]

      idx=bed_Deach$counter %in% paste(bed_Deach$chr,bed_Deach$pos,bed_Deach$ori,sep='/')
      if( !all(idx) ){
        bed_Deach=rbind(bed_Deach,bed_DeachOri[paste(bed_DeachOri$chr,bed_DeachOri$pos,bed_DeachOri$ori,sep='/') %in% bed_Deach$counter[!idx],])
        idx1=bed_Deach$counter %in% paste(bed_Deach$chr,bed_Deach$pos,bed_Deach$ori,sep='/')
        if( !all(idx1)){
          bed_Deach=bed_Deach[bed_Deach$counter %in% paste(bed_Deach$chr,bed_Deach$pos,bed_Deach$ori,sep='/'),]
        }
      }

      if(!all(isEmpty(addtb))){
        bed_Deachadd=addtb[,c(1,2,3,4,which(grepl(paste0('_',m),colnames(addtb))),6)] %>% .[.[,5] >= 1,]
        bed_Deachadd=bed_Deachadd[bed_Deachadd$counter %in% paste(bed_Deach$chr,bed_Deach$pos,bed_Deach$ori,sep='/'),]
      }
    }
    if(grepl('P',AnalysisType) & nrow(bed_P)!=0){
      bed_Peach=bedPmerge[,c(1,2,3,which(grepl(paste0('_',m),colnames(bedPmerge))),6)] %>% .[.[,4] >= 1,]
      bed_Peach=bed_Peach[bed_Peach$counter %in% paste(bed_Peach$chr,bed_Peach$pos,bed_Peach$ori,sep='/'),]
    }
    colnames(bed_Deach)=colnames(bed_Deachadd)=colnames(bed_Peach)=c(c("chr","pos","ori","spl","dis","counter"))

    OutputList[[m]]=rbind(bed_Deach,bed_Deachadd,bed_Peach) %>% unique() %>% .[!duplicated(paste(.$chr,.$pos,.$ori)),]
  }

  if(!any(grepl('chr',chr))){  #update v0.1.6
    for(m in seq_len(length(CaseBam))){
      tmpOut=OutputList[[m]]
      tmpOut$chr=paste0('chr',tmpOut$chr)
      tmpOut$counter=paste0('chr',tmpOut$counter)
      OutputList[[m]]=tmpOut
    }
  }

  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Output generation"))
  for(m in seq_len(length(CaseBam))){
    OutputGenFun(OutputList[[m]],OutputPath,TestID[m])
  }
  writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Done"))
}
