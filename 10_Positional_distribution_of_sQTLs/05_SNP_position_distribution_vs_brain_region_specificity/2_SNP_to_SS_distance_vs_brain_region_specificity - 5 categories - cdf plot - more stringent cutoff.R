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

outputpath="/output/path"
summaryinput="/path/to/summary"
sQTLinput="/path/to/sQTL_run"

library(ggplot2)
library(reshape)

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

for (splicetype in splicetypelist){
  for (type in typelist){
    exoninfopath=paste("/path/to/input_splicing",
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
    
    #2. get a more stringent list of sQTL events
    cutofftype="FDR10"                                                                 #######change########
    if (cutofftype=="FDR1" || cutofftype=="FDR5"){
      setwd(paste("/path/to/result/from/previous/step/",cutofftype,sep=""))
    }
    if (cutofftype=="FDR10"){
      setwd("/path/to/result/from/previous/step/")
    }
    pvmatrix200k=read.table("pvmatrix_200kb.txt",sep="\t",header=T)
    newsQTLexon=rep(NA,dim(pvmatrix200k)[1])
    for (i in 1:length(newsQTLexon)){
      newsQTLexon[i]=strsplit(rownames(pvmatrix200k)[i],split="~")[[1]][1]
    }
    
    #3. get other information based on the new list of sQTL events
    setwd("/path/to/result/from/previous/step/")
    result=read.table(paste(splicetype,"_",type,"_exonclass.txt",sep=""),sep="\t",header=T,check.names=F)
    result=result[newsQTLexon,]
    
    exonclass=as.matrix(result[,1:13])
    num.sig.region=as.numeric(as.matrix(result[,14]))
    
    dis_2_SS=rep(NA,length(newsQTLexon))      #the best category for each sQTL exon
    for (e in 1:length(newsQTLexon)){
      dis_2_SS[e]=pick_category(exonclass[e,],importance_order)[1]
    }
    
    #4. make the plot
    data2plot=cbind(dis_2_SS,num.sig.region)
    df=as.data.frame(data2plot)
    df$dis_2_SS=factor(df$dis_2_SS, levels = importance_order) 
    df$num.sig.region=as.numeric(as.character(df$num.sig.region))

    ###cdf plot version 2###
    setwd(outputpath)
    pdf(paste(splicetype,"_",type,"_SNP_to_SS_distance_vs_brain_region_specificity_CDF_plot2_",cutofftype,".pdf",sep=""),height=6,width=6)
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






