#table of content
#1. stacked bar plot of the number of sQTLs in each brain region (SE+A3SS+A5SS)
#2. stacked bar plot of the number of sQTLs significant in 1, 2, 3, ..., 13 brain regions (SE+A3SS+A5SS)

######################################################################################
#1. stacked bar plot of the number of sQTLs exons in each brain region (SE+A3SS+A5SS)#
######################################################################################
library("data.table")
library(gaston)
library(boot)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
#library(cowplot)

splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
counttype="JC"
PSItype="logit"
outputpath="/output/path"

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

for (type in typelist){
  #read in the result
  totalsqtltable=matrix(0,length(splicetypelist),length(brainregionlist))
  rownames(totalsqtltable)=splicetypelist
  colnames(totalsqtltable)=formalbrainregionlist
  for (splicetype in splicetypelist){
    summaryinput=paste("/path/to/summary",PSItype,counttype,splicetype,sep="/")
    setwd(summaryinput)
    totalsqtlnum=read.table(paste("num_of_sQTL_and_GWASexon_and_diseaseGWASexon_per_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)[2,]
    totalsqtltable[splicetype,]=as.matrix(totalsqtlnum)
  }
  
  #make the bar plot
  df=melt(totalsqtltable)
  colnames(df)=c("Splice.type","Brain.region","Num.of.sQTL.exon")
  df$Splice.type=factor(df$Splice.type, levels=c("A3SS","A5SS","SE"))
  df$Brain.region=factor(df$Brain.region, levels=rev(formalbrainregionlist))
  setwd(outputpath)
  pdf(paste("num_of_sQTL_exon_per_region_",PSItype,"_",counttype,"_",type,".pdf",sep=""))
  if (type=="pvalue"){
    ylimit=900
    breaklist=seq(0,ylimit,by=100)
  }
  if (type=="permutation"){
    ylimit=800
    breaklist=seq(0,ylimit,by=100)
  }
  p=ggplot(data=df, aes(x=Brain.region, y=Num.of.sQTL.exon, fill=Splice.type)) + coord_flip() + scale_y_continuous(limits=c(0,ylimit),breaks=breaklist) + 
    geom_bar(stat="identity") + 
    #scale_fill_brewer(palette="Dark2") + 
    #scale_fill_manual(values=c("#B22222", "#FFA07A", "#FFA500")) + 
    scale_fill_manual(values=c("#BF7DC6", "#E69F00", "#56B4E9")) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      #axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(angle = 60, hjust=1,vjust=1, size=14),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=10),
      axis.line = element_line(colour = "black"),
      
      #background
      panel.background = element_blank(),
      #panel.grid.minor = element_line(colour = "grey"),
      
      #legend
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", colour = "black"),
      legend.text = element_text(size = 10, face = "bold"),
      
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)
  print(p)
  dev.off()
}


################################################################################################################
#2.  stacked bar plot of the number of sQTLs exons significant in 1, 2, 3, ..., 13 brain regions (SE+A3SS+A5SS)#
################################################################################################################
library("data.table")
library(gaston)
library(boot)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
#library(cowplot)

splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
counttype="JC"
PSItype="logit"
outputpath="/output/path"
rootinput=paste("/path/to/sQTL_run",PSItype,counttype,sep="/")
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

for (type in typelist){
  table2plot=matrix(0,length(splicetypelist),length(brainregionlist))
  rownames(table2plot)=splicetypelist
  colnames(table2plot)=as.character(1:length(brainregionlist))
  for (splicetype in splicetypelist){
    sQTLexonlist=c()
    for (b in 1:length(brainregionlist)){
      br=brainregionlist[b]
      inputpath=paste(rootinput,"/",splicetype,"/",br,sep="")
      setwd(inputpath)
      sqtlinfo=read.table(paste("selected.sQTL.glmm.Glimmps_each_exon_cis_",br,".",type,".txt",sep=""),sep="\t")      #sQTL result (exon-SNP pair)
      print(length(unique(as.character(sqtlinfo[,1]))))        #we use sQTL exons
      sQTLexonlist=c(sQTLexonlist,unique(as.character(sqtlinfo[,1])))
    }
    table2plot[splicetype,names(table(table(sQTLexonlist)))]=table(table(sQTLexonlist))
  }
  
  #make the bar plot
  df=melt(table2plot)
  colnames(df)=c("Splice.type","Num.of.brain.region","Num.of.sQTL.exon")
  df$Splice.type=factor(df$Splice.type, levels=c("A3SS","A5SS","SE"))
  df$Num.of.brain.region=factor(df$Num.of.brain.region)
  setwd(outputpath)
  pdf(paste("num_of_sQTL_exon_specificity_",PSItype,"_",counttype,"_",type,".pdf",sep=""))
  if (type=="pvalue"){
    ylimit=1500
    breaklist=seq(0,ylimit,by=100)
  }
  if (type=="permutation"){
    ylimit=900
    breaklist=seq(0,ylimit,by=100)
  }
  p=ggplot(data=df, aes(x=Num.of.brain.region, y=Num.of.sQTL.exon, fill=Splice.type)) + scale_y_continuous(limits=c(0,ylimit),breaks=breaklist) + 
    geom_bar(stat="identity") + 
    #scale_fill_brewer(palette="Dark2") + 
    #scale_fill_manual(values=c("#B22222", "#FFA07A", "#FFA500")) + 
    scale_fill_manual(values=c("#BF7DC6", "#E69F00", "#56B4E9")) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      #axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust=1,vjust=1, size=14),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
      axis.line = element_line(colour = "black"),
      
      #background
      panel.background = element_blank(),
      #panel.grid.minor = element_line(colour = "grey"),
      
      #legend
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", colour = "black"),
      legend.text = element_text(size = 10, face = "bold"),
      
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)
  print(p)
  dev.off()
}



