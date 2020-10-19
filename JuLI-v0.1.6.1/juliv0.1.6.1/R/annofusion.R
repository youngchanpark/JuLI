#' A function for annotation and visualization ver20171204
#'
#' Annotation and visualization.
#' @param Defaults to NULL
#' @keywords annovisual()
#' @export
#' @examples
#' annovisual()


annofusion=function(Output=NULL,
                    Refgene=NULL,
                    Cosmic=NULL,
                    Pfam=NULL,
                    Uniprot=NULL){
  if(isEmpty(Output)){
    MessageFun(5)
  }else{
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Annotation and visualization start"))
    TestID=gsub('^.+/','',Output) %>% gsub('.txt$','',.)
    OutputPath=gsub(paste0('/',gsub('^.+/','',Output)),'',Output)
    dat=fread(Output,sep='\t') %>% setNames(paste0('V',c(1:ncol(.))))
    
    if(nrow(dat)==0){
      MessageFun(3)
    }else{
      if(!isEmpty(Refgene) & !isEmpty(Cosmic) & !isEmpty(Pfam) & !isEmpty(Uniprot)){
        ref= fread(Refgene,showProgress = F) %>% setNames(paste0('V',c(1:ncol(.))))
        cos=read.table(Cosmic,sep='\t',stringsAsFactors = F,header=T)
        pfam=read.table(Pfam,sep='\t',stringsAsFactors = F,fill=T, comment.char = "%")
        gpmap=read.table(Uniprot,sep='\t',stringsAsFactors = F,fill=T,header=T,quote = "")
      }else{
        MessageFun(4)
      }
      
      ff=list()
      
      for(i in 1:nrow(dat)){
        a=reffun(dat$V12[i],dat$V2[i],ref)
        b=reffun(dat$V14[i],dat$V7[i],ref)
        if( dat$V3[i] != dat$V8[i] ){
          if( dat$V13[i] =='+' ){
            gene5p=c(dat$V12[i],dat$V14[i])[c(dat$V3[i],dat$V8[i])==1]
            gene3p=c(dat$V12[i],dat$V14[i])[c(dat$V3[i],dat$V8[i])==0]
            info5p=c(a,b)[c(dat$V3[i],dat$V8[i])==1]
            info3p=c(a,b)[c(dat$V3[i],dat$V8[i])==0]
          }
          if( dat$V13[i] =='-' ){
            gene5p=c(dat$V12[i],dat$V14[i])[c(dat$V3[i],dat$V8[i])==0]
            gene3p=c(dat$V12[i],dat$V14[i])[c(dat$V3[i],dat$V8[i])==1]
            info5p=c(a,b)[c(dat$V3[i],dat$V8[i])==0]
            info3p=c(a,b)[c(dat$V3[i],dat$V8[i])==1]
          }
        }
        if( dat$V3[i] == dat$V8[i] ){
          if( dat$V3[i] ==1 ){
            gene5p=c(dat$V12[i],dat$V14[i])[c(dat$V13[i],dat$V15[i])=='+']
            gene3p=c(dat$V12[i],dat$V14[i])[c(dat$V13[i],dat$V15[i])=='-']
            info5p=c(a,b)[c(dat$V13[i],dat$V15[i])=='+']
            info3p=c(a,b)[c(dat$V13[i],dat$V15[i])=='-']
          }
          if( dat$V3[i] ==0 ){
            gene5p=c(dat$V12[i],dat$V14[i])[c(dat$V13[i],dat$V15[i])=='-']
            gene3p=c(dat$V12[i],dat$V14[i])[c(dat$V13[i],dat$V15[i])=='+']
            info5p=c(a,b)[c(dat$V13[i],dat$V15[i])=='-']
            info3p=c(a,b)[c(dat$V13[i],dat$V15[i])=='+']
          }
        }
        c=paste0(gene5p,'->',gene3p)
        d=NA
        if(grepl('Intron',a) & grepl('Intron',b)){
          if( gsub('^.+Frame','',a)==gsub('^.+Frame','',b) ){d='Inframe'}else{d='Outframe'}
        }
        if(sum(grepl('Exon',c(a,b)))==1){
          if( grepl('Intron',info5p) ){
            if( gsub('^.+Frame','',info5p)==gsub('[[:digit:]]+,','',gsub('^.+Frame','',info3p)) ){d='Possible_Inframe'}else{d='Possible_Outframe'}
          }
          if( grepl('Exon',info5p) ){d='Possible_Outframe'}
        }
        e=NA
        if( dat$V12[i]!=dat$V14[i] ){
          if( sum(grepl(dat$V12[i],cos$Translocation.Name) & grepl(dat$V14[i],cos$Translocation.Name))!=0 ){
            g=cos$Translocation.Name[grepl(dat$V12[i],cos$Translocation.Name) & grepl(dat$V14[i],cos$Translocation.Name)]
            g1=names(sort(table(gsub('[{].+$','',g)),decreasing=T))[1]
            if( gsub('->.+$','',c)==g1 ){
              e1=sort(table(cos$Primary.site[grepl(dat$V12[i],cos$Translocation.Name) & grepl(dat$V14[i],cos$Translocation.Name)]),decreasing=T)
              e=paste(paste0(names(e1),'(',e1,')'),collapse ='_')
            }
          }
        }
        ff[[i]]=data.frame(a,b,c,d,e,stringsAsFactors = F)
      }
      
      t.dat=data.frame(dat,rbindlist(ff),stringsAsFactors = F)
      colnames(t.dat)=c('ChrA','BreakA','OriA','DisA','SplitA','ChrB','BreakB','OriB','DisB','SplitB','Event','GeneA','StrGeneA','GeneB','StrGeneB','InfoA','InfoB','Direction','Frame','Cosmic')
      
      tt.dat=t.dat[order(t.dat$DisA+t.dat$DisB,t.dat$SplitA+t.dat$SplitB,decreasing=T),]
      ft.dat=tt.dat[ !grepl('Compl',tt.dat$GeneA) & !grepl('Compl',tt.dat$GeneB) & !grepl('Flanking',tt.dat$GeneA) & !grepl('Flanking',tt.dat$GeneB) & grepl('NM',tt.dat$InfoA) & grepl('NM',tt.dat$InfoB),] #gene-gene fusion
      
      fwrite(tt.dat,paste0(OutputPath,'/',TestID,'.annotated.txt'),sep='\t',showProgress=F)
#     fwrite(ft.dat,paste0(OutputPath,'/',TestID,'.annotated.gene.vcf'),sep='\t',showProgress=F)
      
      if(nrow(ft.dat)!=0 ){
        plots = list()
        for(p in 1:nrow(ft.dat)){
          g5=unlist(strsplit(ft.dat$Direction[p],'->'))[1]
          g3=unlist(strsplit(ft.dat$Direction[p],'->'))[2]
          if(g5!=g3){
            info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g5][1]
            info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g3][1]
            bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g5][1]
            bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$GeneA[p],ft.dat$GeneB[p])==g3][1]
          }
          if(g5==g3){
            if(ft.dat$StrGeneA[p]=='+'){
              info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
            }
            if(ft.dat$StrGeneA[p]=='-'){
              info5=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              info3=c(ft.dat$InfoA[p],ft.dat$InfoB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              bp5=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]
              bp3=c(ft.dat$BreakA[p],ft.dat$BreakB[p])[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))]        
            }
          }
          
          gd5=data.frame(0.4,mapfun(g5,info5,bp5,gpmap,pfam,ref),stringsAsFactors=F)
          gd3=data.frame(0.4,mapfun(g3,info3,bp3,gpmap,pfam,ref),stringsAsFactors=F)
          gd3$st=gd3$st+0.5
          gd3$en=gd3$en+0.5
          gd3$nbp=gd3$nbp+0.5
          colnames(gd5)[1]=colnames(gd3)[1]='V1'
          
          fgd5=gd5[c(1:as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))),]
          if(!grepl('Intron',info5)){
            fgd5$en[as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))]=fgd5$nbp[1]
            no5=0.5-fgd5$en[nrow(fgd5)]
            fgd5$st=fgd5$st+no5
            fgd5$en=fgd5$en+no5
          }
          if(grepl('Intron',info5)){
            no5=(0.5-(fgd5$nbp[1]-fgd5$en[as.numeric(gsub('^.+[(]','',gsub('/.+$','',info5)))]))-fgd5$en[nrow(fgd5)]
            fgd5$st=fgd5$st+no5
            fgd5$en=fgd5$en+no5
          }
          fgd3=gd3[c(as.numeric(gsub('^.+[(]','',gsub('/.+$','',info3))):nrow(gd3)),]
          if(!grepl('Intron',info3)){
            fgd3$st[1]=fgd3$nbp[1]
            no3=fgd3$st[1]-0.5
            fgd3$st=fgd3$st-no3
            fgd3$en=fgd3$en-no3
          }
          if(grepl('Intron',info3)){
            fgd3=fgd3[-1,]
            no3=fgd3$nbp[1]-0.5
            fgd3$st=fgd3$st-no3
            fgd3$en=fgd3$en-no3
          }      
          fgd=rbind(fgd5,fgd3)
          fgd$V1=0.2
          
          gd=rbind(gd5,gd3,fgd)  
          if(sum(!is.na(gd$Domain))==0){gd$Domain='NA'}
          gd$Domain=factor(gd$Domain)
          
          arr1=c(fgd5$nbp[1],0.498,0.352,0.25)
          arr2=c(fgd3$nbp[1],0.502,0.352,0.25)
          arr.dat=data.frame(rbind(arr1,arr2),stringsAsFactors = F)
          
          text0=c(0.05,1,paste('ID:',TestID))     
          text1=c(0.05,0.95,paste('Fusion:',ft.dat$Direction[p]))
          text2=c(0.05,0.9,paste('Frame:',ft.dat$Frame[p]))
          text3=c(0.05,0.85,paste('Event type:',ft.dat$Event[p]))
          text4=c(0.05,0.8,paste('Supporting reads:',paste0('split=',ft.dat$SplitA[p]+ft.dat$SplitB[p],'/','discordant=',ft.dat$DisA[p]+ft.dat$DisB[p])))
          text5=c(0.05,0.75,paste('Cosmic:',ft.dat$Cosmic[p]))
          text6=c(0.05,0.7,'GeneName/Breakpoint/Strand/FusionInfo:')
          InA=paste0(' ',ft.dat$GeneA[p],'/',ft.dat$ChrA[p],':',ft.dat$BreakA[p],'/',paste0('(',ft.dat$StrGeneA[p],')'),'/',ft.dat$InfoA[p])
          InB=paste0(' ',ft.dat$GeneB[p],'/',ft.dat$ChrB[p],':',ft.dat$BreakB[p],'/',paste0('(',ft.dat$StrGeneB[p],')'),'/',ft.dat$InfoB[p])
          if(g5!=g3){
            text7=c(0.05,0.65,c(InA,InB)[c(ft.dat$InfoA[p],ft.dat$InfoB[p])==info5][1])
            text8=c(0.05,0.6,c(InA,InB)[c(ft.dat$InfoA[p],ft.dat$InfoB[p])==info3][1])
          }
          if(g5==g3){
            if(ft.dat$StrGeneA[p]=='+'){
              text7=c(0.05,0.65,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
              text8=c(0.05,0.6,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
            }
            if(ft.dat$StrGeneA[p]=='-'){
              text7=c(0.05,0.65,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==max(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
              text8=c(0.05,0.6,c(InA,InB)[c(ft.dat$BreakA[p],ft.dat$BreakB[p])==min(c(ft.dat$BreakA[p],ft.dat$BreakB[p]))])
            }
          }
          text9=c(0.05,0.05,paste0(p,'/',nrow(ft.dat)))
          text10=c(0.05,0.46,paste('5`',g5))
          text11=c(0.55,0.46,paste('5`',g3))
          text12=c(0.45,0.12,ft.dat$Direction[p])
          
          text.dat=data.frame(rbind(text0,text1,text2,text3,text4,text5,text6,text7,text8,text9,text10,text11,text12),stringsAsFactors = F)
          text.dat$X1=as.numeric(text.dat$X1)
          text.dat$X2=as.numeric(text.dat$X2)
          
          plots[[p]]=ggplot() +
            geom_segment(data=gd,aes(x=st, xend=en, y=V1, yend=V1,colour=Domain), size=15) +
            scale_colour_hue(l=40) +
            geom_segment(data=arr.dat,aes(x = X1, y = X3, xend = X2, yend = X4),arrow = arrow(length = unit(0.02, "npc"))) +
            geom_text(data=text.dat,aes(x=X1, y=X2, label=X3),hjust=0)+
            ylim(c(0,1)) +
            theme(axis.ticks=element_blank(),axis.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
                  legend.position="right",legend.key = element_blank(),legend.background =element_blank()) 
        }
        
        pdf(paste0(OutputPath,'/',TestID,'.annotated.gene.pdf'),width=10,height=8)
        suppressWarnings(invisible(lapply(plots, print)))
        invisible(dev.off())
      }
    }
    writeLines(paste0("[",format(Sys.time(),"%Y-%b-%d %H:%M:%S"),"] ","Annotation and visualization finish"))
  }
}
