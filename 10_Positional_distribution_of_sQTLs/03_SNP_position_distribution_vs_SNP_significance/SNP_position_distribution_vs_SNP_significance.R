library(ggplot2)
library(ggsci)
library(scales)
splicetypelist=c("SE","A3SS","A5SS")                
typelist=c("pvalue","permutation")

setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/scripts")
brainregion=as.character(as.matrix(read.table("brainregionlist.txt",sep="\t")))
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC/"
rootsQTLinput='/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC'
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/SNP_position_distribution"

findpos=function(x,exoninfo,searchrange,splicetype){
  strand=strsplit(exoninfo,split="\\|")[[1]][5]
  SNPpos=as.numeric(strsplit(as.matrix(x["SNPID"]),split="_")[[1]][2])
  
  if (splicetype=="SE"){
    exonstart=as.numeric(strsplit(exoninfo,split="\\|")[[1]][6])
    exonend=as.numeric(strsplit(exoninfo,split="\\|")[[1]][7])
    if (strand=="+"){
      SS5start=exonend-3
      SS5end=exonend+6
      SS3start=exonstart-20
      SS3end=exonstart+3
    }
    if (strand=="-"){
      SS5start=exonstart-6
      SS5end=exonstart+3
      SS3start=exonend-3
      SS3end=exonend+20
    }

    dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end)) 
    
    if (SNPpos>=SS5start && SNPpos<=SS5end){ 
      label="5'SS"
    }else if (SNPpos>=SS3start && SNPpos<=SS3end){
      label="3'SS"
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
    }else if (dis<=searchrange){
      label="<=300bp"
    }else {
      label=">300bp"
    }
  }

  if (splicetype=="A5SS"){
    exonstart=as.numeric(strsplit(exoninfo,split="\\|")[[1]][6])       #we choose the longer exon
    exonend=as.numeric(strsplit(exoninfo,split="\\|")[[1]][7])
    if (strand=="+"){
      ass1=as.numeric(strsplit(exoninfo,split="\\|")[[1]][7])
      ass1_start=ass1-3
      ass1_end=ass1+6
      ass2=as.numeric(strsplit(exoninfo,split="\\|")[[1]][9])
      ass2_start=ass2-3
      ass2_end=ass2+6
    }
    if (strand=="-"){
      ass1=as.numeric(strsplit(exoninfo,split="\\|")[[1]][6])
      ass1_start=ass1-6
      ass1_end=ass1+3
      ass2=as.numeric(strsplit(exoninfo,split="\\|")[[1]][8])
      ass2_start=ass2-6
      ass2_end=ass2+3
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="5'SS"
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="5'SS"
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
    }else if (dis<=searchrange){
      label="<=300bp"
    }else {
      label=">300bp"
    }
  }
  
  if (splicetype=="A3SS"){
    exonstart=as.numeric(strsplit(exoninfo,split="\\|")[[1]][6])       #we choose the longer exon
    exonend=as.numeric(strsplit(exoninfo,split="\\|")[[1]][7])
    if (strand=="+"){
      ass1=as.numeric(strsplit(exoninfo,split="\\|")[[1]][6])
      ass1_start=ass1-20
      ass1_end=ass1+3
      ass2=as.numeric(strsplit(exoninfo,split="\\|")[[1]][8])
      ass2_start=ass2-20
      ass2_end=ass2+3
    }
    if (strand=="-"){
      ass1=as.numeric(strsplit(exoninfo,split="\\|")[[1]][7])
      ass1_start=ass1-3
      ass1_end=ass1+20
      ass2=as.numeric(strsplit(exoninfo,split="\\|")[[1]][9])
      ass2_start=ass2-3
      ass2_end=ass2+20
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="3'SS"
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="3'SS"
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
    }else if (dis<=searchrange){
      label="<=300bp"
    }else {
      label=">300bp"
    }
  }  
  
  return(label)
}


for (type in typelist){
  for (splicetype in splicetypelist){
    #the matrix which stores the result of all brain regions
    allposmatrix=matrix(,nrow=0,ncol=4)
    colnames(allposmatrix)=c("sQTL.exon","SNP","p.value","category")
    
    for (i in 1:length(brainregion)[1]){
      br=brainregion[i]
      print(br)
      inputpath=paste(rootinput,splicetype,sep="")
      #read in the sQTL result in the current brain region
      setwd(inputpath)
      exonSNP=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",br,"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)    #the results based on pvalue and permutation are the same
      
      if (dim(exonSNP)[1]>0){
        #get list of sQTL exons
        sQTLexon=unique(exonSNP[,c("exon_full_ID","exon_short_ID")])
        
        #for each exon, get all the SNPs within 200kb of that exon and group each SNP into the 5 categories
        positionmatrix=matrix(,nrow=0,ncol=4)
        colnames(positionmatrix)=c("sQTL.exon","SNP","p.value","category")
        
        for (e in 1:dim(sQTLexon)[1]){
          sQTLexonpath=paste(rootsQTLinput,splicetype,br,paste("Glimmps_each_exon_cis_",br,sep=""),sep="/")
          setwd(sQTLexonpath)
          singleexonresult=read.table(paste(as.character(sQTLexon[e,"exon_short_ID"]),".asso",sep=""),sep="\t",header=T)
          pos=apply(singleexonresult,1,findpos,exoninfo=as.character(sQTLexon[e,"exon_full_ID"]),searchrange=300, splicetype)
          temp=cbind(rep(as.character(sQTLexon[e,"exon_full_ID"]),dim(singleexonresult)[1]),as.matrix(singleexonresult[,"SNPID"]),-log10(singleexonresult[,"pvals.lm"]),pos)
          positionmatrix=rbind(positionmatrix,temp)
        }
        
        #make the box plot
        df=as.data.frame(positionmatrix[,c("p.value","category")])
        df$p.value=as.numeric(as.matrix(df$p.value))
        df$category <- factor(df$category, levels = c("5'SS","3'SS","exon","<=300bp",">300bp"))   #fix the order of diseases
        
        setwd(outputpath)
        pdf(paste(splicetype,type,br,"SNP position distribution.pdf"),height=6,width=6)
        p <- ggplot(df, aes(x=category, y=p.value, fill=category)) + 
          geom_boxplot() + 
          labs(y = expression('-log'[10]*'(P)'), x="") + 
          ggtitle(br) + 
          scale_fill_manual(values=c("#913943", "#0E730C", "#0F0C73","#26908E","#626565")) +
          theme(
            # axis
            axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
            axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
            axis.text = element_text(size = rel(1.1)),
            axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
            axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
            axis.line = element_line(colour = "black"),
            
            #background
            panel.background = element_blank(),
            
            #legend
            legend.position="none",
            
            # strip
            strip.text=element_text(size = rel(1.3)),
            aspect.ratio=1,
            complete = T)
        print(p)
        dev.off()
        
        #add the result in the current brain region to all result
        allposmatrix=rbind(allposmatrix,positionmatrix)
      }
    }
    #generate the same plot using the result from all brain regions
    #make the box plot
    allposmatrix[which(allposmatrix[,"category"] %in% ">300bp"),"category"]=paste("Intron","(>300bp)",sep="\n")
    allposmatrix[which(allposmatrix[,"category"] %in% "<=300bp"),"category"]=paste("Intron","(<=300bp)",sep="\n")
    allposmatrix[which(allposmatrix[,"category"] %in% "exon"),"category"]="Exon"
    df.all=as.data.frame(allposmatrix[,c("p.value","category")])
    df.all$p.value=as.numeric(as.matrix(df.all$p.value))
    df.all$category <- factor(df.all$category, levels = c("5'SS","3'SS","Exon",paste("Intron","(<=300bp)",sep="\n"),paste("Intron","(>300bp)",sep="\n")))   #fix the order of diseases
    
    setwd(outputpath)
    pdf(paste(splicetype,type,"all region SNP position distribution.pdf"),height=6,width=6)
    p <- ggplot(df.all, aes(x=category, y=p.value, fill=category)) + 
      geom_boxplot() + 
      labs(y = expression('-log'[10]*'(P)'), x="") + 
      scale_fill_manual(values=c("#913943", "#0E730C", "#0F0C73","#26908E","#626565")) +
      theme(
        # axis
        axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
        axis.text = element_text(size = rel(1.1)),
        axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
        axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
        axis.line = element_line(colour = "black"),
        
        #background
        panel.background = element_blank(),
        
        #legend
        legend.position="none",
        
        # strip
        strip.text=element_text(size = rel(1.3)),
        aspect.ratio=1,
        complete = T)
    print(p)
    dev.off()
  }
}


