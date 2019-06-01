brainregionlist=c("Brain-Amygdala",
                  "Brain-AnteriorcingulatecortexBA24",
                  "Brain-Caudatebasalganglia",
                  "Brain-CerebellarHemisphere",
                  "Brain-Cerebellum",
                  "Brain-Cortex",
                  "Brain-FrontalCortexBA9",
                  "Brain-Hippocampus",
                  "Brain-Hypothalamus",
                  "Brain-Nucleusaccumbensbasalganglia",
                  "Brain-Putamenbasalganglia",
                  "Brain-Spinalcordcervicalc-1",
                  "Brain-Substantianigra")
formalbrainregionlist=c("Amygdala",
                        "Anterior cingulate cortex BA24",
                        "Caudate basal ganglia",
                        "Cerebellar Hemisphere",
                        "Cerebellum",
                        "Cortex",
                        "Frontal Cortex BA9",
                        "Hippocampus",
                        "Hypothalamus",
                        "Nucleus accumbens basal ganglia",
                        "Putamen basal ganglia",
                        "Spinal cord cervical c-1",
                        "Substantia nigra")

splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
PSItype="logit"
counttype="JC"

outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/SNP_to_SS_distance_vs_brain_region_specificity"
summaryinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary"
sQTLinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"

library(ggplot2)
library(reshape)


findpos=function(exoninfo,shortexonID,snppos,splicetype){         #classify SNPs into five categories
  oneexoninfo=as.character(as.matrix(exoninfo[shortexonID,]))
  strand=oneexoninfo[5]
  exonstart=as.numeric(oneexoninfo[6])
  exonend=as.numeric(oneexoninfo[7])
  SNPpos=as.numeric(snppos)
  
  if (splicetype=="SE"){
    if (strand=="+"){
      #wide splice site
      SS5start=exonend-3
      SS5end=exonend+6
      SS3start=exonstart-20
      SS3end=exonstart+3
      #narrow splice site
      diSS5start=exonend+1
      diSS5end=exonend+2
      diSS3start=exonstart-1
      diSS3end=exonstart
    }
    if (strand=="-"){
      #wide splice site
      SS5start=exonstart-6
      SS5end=exonstart+3
      SS3start=exonend-3
      SS3end=exonend+20
      #narrow splice site
      diSS5start=exonstart-1
      diSS5end=exonstart
      diSS3start=exonend+1
      diSS3end=exonend+2
    }
    newexonstart=exonstart
    newexonend=exonend
    
    dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end)) 
    
    if (SNPpos>=diSS5start && SNPpos<=diSS5end){
      label="dinucleotide"    #5'SS dinucleotide
      return(label)
    }else if (SNPpos>=diSS3start && SNPpos<=diSS3end){
      label="dinucleotide"    #3'SS dinucleotide
      return(label)
    }else if (SNPpos>=SS5start && SNPpos<=SS5end){
      label="SS"      #5'SS
      return(label)
    }else if (SNPpos>=SS3start && SNPpos<=SS3end){
      label="SS"      #3'SS
      return(label)
    }else if (SNPpos>=newexonstart && SNPpos<=newexonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=300bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
  
  if (splicetype=="A3SS"){
    if (strand=="+"){
      ass1=as.numeric(oneexoninfo[6])
      ass1_start=ass1-20
      ass1_end=ass1+3
      ass2=as.numeric(oneexoninfo[8])
      ass2_start=ass2-20
      ass2_end=ass2+3
      #narrow splice site
      diass1start=ass1-1
      diass1end=ass1
      diass2start=ass2-1
      diass2end=ass2
    }
    if (strand=="-"){
      ass1=as.numeric(oneexoninfo[7])
      ass1_start=ass1-3
      ass1_end=ass1+20
      ass2=as.numeric(oneexoninfo[9])
      ass2_start=ass2-3
      ass2_end=ass2+20
      #narrow splice site
      diass1start=ass1+1
      diass1end=ass1+2
      diass2start=ass2+1
      diass2end=ass2+2
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=diass1start && SNPpos<=diass1end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=diass2start && SNPpos<=diass2end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="SS"
      return(label)
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="SS"
      return(label)
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=300bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
  
  if (splicetype=="A5SS"){
    if (strand=="+"){
      ass1=as.numeric(oneexoninfo[7])
      ass1_start=ass1-3
      ass1_end=ass1+6
      ass2=as.numeric(oneexoninfo[9])
      ass2_start=ass2-3
      ass2_end=ass2+6
      #narrow splice site
      diass1start=ass1+1
      diass1end=ass1+2
      diass2start=ass2+1
      diass2end=ass2+2
    }
    if (strand=="-"){
      ass1=as.numeric(oneexoninfo[6])
      ass1_start=ass1-6
      ass1_end=ass1+3
      ass2=as.numeric(oneexoninfo[8])
      ass2_start=ass2-6
      ass2_end=ass2+3
      #narrow splice site
      diass1start=ass1-1
      diass1end=ass1
      diass2start=ass2-1
      diass2end=ass2
    }
    
    dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    
    if (SNPpos>=diass1start && SNPpos<=diass1end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=diass2start && SNPpos<=diass2end){
      label="dinucleotide"
      return(label)
    }else if (SNPpos>=ass1_start && SNPpos<=ass1_end){ 
      label="SS"
      return(label)
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="SS"
      return(label)
    }else if (SNPpos>=exonstart && SNPpos<=exonend){
      label="exon"
      return(label)
    }else if (dis<=300){
      label="<=100bp"
      return(label)
    }else{
      label=">300bp"
      return(label)
    }
  }
}

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

for (splicetype in splicetypelist){
  for (type in typelist){
    exoninfopath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing",
                       PSItype,counttype,splicetype,sep="/")
    setwd(exoninfopath)
    exoninfo=read.table(paste("exon_info.fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)
    rownames(exoninfo)=exoninfo[,"ID"]
    exonshorterID=rownames(exoninfo)
    
    #1. get all sQTL exons
    sQTLexon=c()
    inputpath=paste(summaryinput,PSItype,counttype,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
    }
    sQTLexon=unique(sQTLexon)
    
    #2. for each sQTL exon:
    #(1) find the category of all the significant SNPs and pick the best one in each brain region
    #(2) count the number of significant regions
    num.sig.region=rep(0,length(sQTLexon))
    exonclass=matrix(NA,length(sQTLexon),length(brainregionlist))     #store the category
    colnames(exonclass)=brainregionlist
    rownames(exonclass)=sQTLexon
    sigexoninfo=matrix(0,length(sQTLexon),length(brainregionlist))     #1: significant. 0: insignificant
    colnames(sigexoninfo)=brainregionlist
    rownames(sigexoninfo)=sQTLexon

    for (e in 1:length(sQTLexon)){
      for (br in 1:length(brainregionlist)){
        #get the p value cutoff
        if (type=="pvalue"){
          cutoff=10^-5
        }
        if (type=="permutation"){
          setwd(paste(sQTLinput,"/",PSItype,"/",counttype,"/",splicetype,"/",brainregionlist[br],sep=""))
          cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
        }
        
        shortexonID=paste("SE_",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="")
        path=paste(sQTLinput,PSItype,counttype,splicetype,brainregionlist[br],paste("Glimmps_each_exon_cis_",brainregionlist[br],sep=""),sep="/")
        setwd(path)
        temp=read.table(paste(shortexonID,".asso",sep=""),header=T)
        #get the significant SNP
        subtemp=subset(temp,temp[,"pvals.lm"]<=cutoff)
        if (dim(subtemp)[1]>0){        #if we have significant SNPs
          num.sig.region[e]=num.sig.region[e]+1
          sigexoninfo[e,br]=1
          
          #for all the significant SNPs, we classify them into categories
          subclasslist=rep(NA,dim(subtemp)[1])
          for (j in 1:dim(subtemp)[1]){
            subclasslist[j]=findpos(exoninfo,shortexonID,subtemp[j,"Pos"],splicetype)
          }
          #we find the best category
          exonclass[e,br]=pick_category(subclasslist,importance_order)[1]
        }
      }
    }
    #we may have NA in the matrix, this is because for some exons in some brain regions, there may be no significant SNPs
    setwd(outputpath)
    write.table(cbind(exonclass,num.sig.region),paste(splicetype,"_",type,"_exonclass.txt",sep=""),sep="\t")
    write.table(sigexoninfo,paste(splicetype,"_",type,"_sigexoninfo.txt",sep=""),sep="\t")
    
    dis_2_SS=rep(NA,length(sQTLexon))      #the best category for each sQTL exon
    for (e in 1:length(sQTLexon)){
      dis_2_SS[e]=pick_category(exonclass[e,],importance_order)[1]
    }
    
    #4. make the plot
    data2plot=cbind(dis_2_SS,num.sig.region)
    df=as.data.frame(data2plot)
    df$dis_2_SS=factor(df$dis_2_SS, levels = importance_order) 
    df$num.sig.region=as.numeric(as.character(df$num.sig.region))
    
    ###violin plot###
    setwd(outputpath)
    pdf(paste(splicetype,type,"SNP_to_SS_distance_vs_brain_region_specificity_violin_plot.pdf"),height=6,width=6)
    p <- ggplot(df, aes(x=dis_2_SS, y=num.sig.region, fill=dis_2_SS)) + 
      scale_fill_manual(values=c("#FDC086", "#BEAED4", "#0F0C73","#26908E","#626565")) +
      geom_violin(scale="width") + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + geom_boxplot(width=0.2) + 
      labs(x = "Location of SNP", y="No. of significant regions") + 
      scale_y_continuous(breaks = round(seq(0, 13, by = 1),1)) + 
      theme(
        # axis
        axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
        axis.text = element_text(size = rel(1.1)),
        axis.text.x = element_text(hjust = 0.5, vjust = 0, size=12),
        axis.text.y = element_text(vjust = 0.5, hjust = 0, size=12),
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
    

    ###cdf plot###
    cdfdf=melt(df)
    setwd(outputpath)
    pdf(paste(splicetype,type,"SNP_to_SS_distance_vs_brain_region_specificity_CDF_plot.pdf"),height=4,width=4)
    p=ggplot(cdfdf, aes(x=value)) + stat_ecdf(aes(colour=dis_2_SS)) + 
      #scale_color_brewer(palette="Accent") +
      scale_color_manual(values=c("#FDC086", "#BEAED4", "#0F0C73","#26908E","#626565")) + 
      #scale_color_npg() + 
      labs(x = "number of significant regions", y="") +
      scale_x_continuous(breaks = round(seq(0, 13, by = 1),1)) +  
      theme(
        # axis
        axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
        axis.text = element_text(size = rel(1.1)),
        axis.text.x = element_text(hjust = 0.5, vjust = 0, size=8),
        axis.text.y = element_text(vjust = 0.5, hjust = 0, size=12),
        axis.line = element_line(colour = "black"),
        
        #legend
        legend.key=element_blank(),
        
        #background
        panel.background = element_blank(),
        
        # strip
        strip.text=element_text(size = rel(1.3)),
        aspect.ratio=1,
        complete = T)
    print(p)
    dev.off()
    
    
    ###cdf plot version 2###
    setwd(outputpath)
    pdf(paste(splicetype,type,"SNP_to_SS_distance_vs_brain_region_specificity_CDF_plot2.pdf"),height=6,width=6)
    plot(ecdf(df[df$dis_2_SS=="dinucleotide",]$num.sig.region),
         verticals=TRUE, 
         do.points=FALSE,
         xlim=c(1,13),
         xlab="Number of significant regions",
         ylab="Percent",
         xaxt = 'n',
         main = "",
         lwd=5,
         bty="n",
         #xaxs="i",
         #yaxs="i",
         col="#FDC086")
    axis(1, at=1:13, labels=1:13)
    lines(ecdf(df[df$dis_2_SS=="SS",]$num.sig.region),
          verticals=TRUE, do.points=FALSE,col.01line = NULL,
          lwd=5,
          xlim=c(1,13),
          col="#BEAED4")
    lines(ecdf(df[df$dis_2_SS=="exon",]$num.sig.region),
          verticals=TRUE, do.points=FALSE,col.01line = NULL,
          lwd=5,
          xlim=c(1,13),
          col="#0F0C73")
    lines(ecdf(df[df$dis_2_SS=="<=300bp",]$num.sig.region),
          verticals=TRUE, do.points=FALSE,col.01line = NULL,
          lwd=5,
          xlim=c(1,13),
          col="#26908E")
    lines(ecdf(df[df$dis_2_SS==">300bp",]$num.sig.region),
          verticals=TRUE, do.points=FALSE,col.01line = NULL,
          lwd=5,
          xlim=c(1,13),
          col="#626565")
    dev.off()
  }
}







