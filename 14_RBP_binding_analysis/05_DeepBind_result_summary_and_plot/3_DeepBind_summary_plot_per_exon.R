#purpose of this plot:
#for each exon, we get its deepbind score (DBscore + WT and MUT score) for all the snps and all the RBPs
#then we make the bubble plot using these values

args <- commandArgs(TRUE)
splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
type="pvalue"       #pvalue or permutation
shortID=args[3]
#shortID="SE_189436"    #MAPT
#shortID="SE_128125"    #CXXC5

windowsize=20

library(NMF)
library(RColorBrewer)
require(ggplot2)
require(ggseqlogo)
#require(cowplot)
#library(reshape)
require(scales)
library(reshape2)
library(gridExtra)
library(grid)
#library(ComplexHeatmap)
#library(circlize)
#library(lemon)
library(ggsci)
require("ggrepel")

#get all sQTL exons
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC"
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
sQTLexon=c()
inputpath=paste(rootinput,splicetype,sep="/")
setwd(inputpath)
for (i in 1:length(brainregionlist)){
  temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
  sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
}
sQTLexon=unique(sQTLexon)
shortIDlist=rep(NA,length(sQTLexon))
for (e in 1:length(sQTLexon)){
  shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
}

#read in the RBP information
RBPdb="/u/nobackup/yxing/PROJECT/yidazhan/research/software/deepbind/db/db.tsv"
RBPtable=read.table(RBPdb,sep="\t",header=T)
subRBPtable=subset(RBPtable,RBPtable[,"Species"]=="Homo sapiens")
subRBPtable=subset(subRBPtable,subRBPtable[,"Type"]=="RBP")

#read in the joblist
rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(rootoutput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")

#generate output folder
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_result_plot_per_exon/",splicetype,"/",type,sep="")
command=paste("mkdir -p",outputpath)
system(command)

#get the information for all the jobs related to this exon
exonjoblist=subset(uniquejoblist,uniquejoblist[,"Exon"] %in% shortID)

#root folder of the current exon
fullexonID=gsub("\\|",",",sQTLexon[which(shortIDlist %in% shortID)])
rootexonresult=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_output/",
                     splicetype,"/",type,"/",fullexonID,sep="")

###################################################################
#read in the DB score for all the jobs related to the current exon#
###################################################################
DBscore_WT=DBscore_MUT=DBscore=rep(0,dim(exonjoblist)[1])
data2plot=data2plot0=data2plot00=cbind(exonjoblist,DBscore_WT,DBscore_MUT,DBscore)      #three different versions of score

for (i in 1:dim(exonjoblist)[1]){
  print(i)
  subfolder=paste(exonjoblist[i,"SNP"],exonjoblist[i,"RBP"],sep="~")
  setwd(paste(rootexonresult,subfolder,sep="/"))
  
  #version 1: original score
  dbscore=read.table("dbscore.txt",header=T,sep="\t",check.names=F)
  dbscore_WT=read.table("dbscore_WT.txt",header=T,sep="\t",check.names=F)
  dbscore_MUT=read.table("dbscore_MUT.txt",header=T,sep="\t",check.names=F)
  
  #get reference and alternative allele
  snpid=as.character(exonjoblist[i,"SNP"])
  allele0=strsplit(snpid,split="_")[[1]][3]     #reference allele, the value is 0
  allele1=strsplit(snpid,split="_")[[1]][4]     #alternative allele, the value is 1
  if (allele0==colnames(dbscore)[windowsize+1]){          #the SNP is on the same strand with the exon sequence
    ref=allele0
    alt=allele1
  }else{                                        #the SNP is on the opposite strand of the exon sequence
    ref=chartr("ATGC","TACG",allele0)           #reverse complement
    alt=chartr("ATGC","TACG",allele1)
  }
  
  data2plot[i,"DBscore"]=dbscore[alt,c(windowsize+1)]  
  data2plot[i,"DBscore_WT"]=dbscore_WT[alt,c(windowsize+1)]  
  data2plot[i,"DBscore_MUT"]=dbscore_MUT[alt,c(windowsize+1)]  
  
  
  #version 2: at least one positive score
  dbscore=read.table("dbscore.txt",header=T,sep="\t",check.names=F)
  dbscore_WT0=read.table("dbscore_WT0.txt",header=T,sep="\t",check.names=F)
  dbscore_MUT0=read.table("dbscore_MUT0.txt",header=T,sep="\t",check.names=F)
  
  data2plot0[i,"DBscore"]=dbscore[alt,c(windowsize+1)]  
  data2plot0[i,"DBscore_WT"]=dbscore_WT0[alt,c(windowsize+1)]  
  data2plot0[i,"DBscore_MUT"]=dbscore_MUT0[alt,c(windowsize+1)]  
  
  
  #version 3: both positive score
  dbscore=read.table("dbscore.txt",header=T,sep="\t",check.names=F)
  dbscore_WT00=read.table("dbscore_WT00.txt",header=T,sep="\t",check.names=F)
  dbscore_MUT00=read.table("dbscore_MUT00.txt",header=T,sep="\t",check.names=F)
  
  data2plot00[i,"DBscore"]=dbscore[alt,c(windowsize+1)]  
  data2plot00[i,"DBscore_WT"]=dbscore_WT00[alt,c(windowsize+1)]  
  data2plot00[i,"DBscore_MUT"]=dbscore_MUT00[alt,c(windowsize+1)]  
}

##########################
#generate the bubble plot#
##########################
bubbleplot=function(data,path,fullexon,version){
  hi=6
  wi=6
  upperlim=ceiling(max(range(data[,"DBscore_WT"]),range(data[,"DBscore_MUT"])))
  lowerlim=floor(min(range(data[,"DBscore_WT"]),range(data[,"DBscore_MUT"])))
  setwd(path)
  pdf(paste(fullexonID,"_withtext_",version,".pdf",sep=""),height=hi,width=wi)
  p=ggplot(data, aes(x=DBscore_WT, y=DBscore_MUT, size=abs(DBscore), color=SNP, label=RBP)) +
    geom_point(alpha=0.7) +
    #scale_size_continuous( trans="exp", range=c(1, 18)) +     # highlight very high variables (exponential transformation of original value)
    #geom_label_repel(aes(label = RBP,
    #                     fill = factor(SNP)), color = 'white',
    #                 size = 1)  + 
    geom_text(size=1,aes(color=factor(SNP))) + 
    #scale_color_npg() + 
    scale_color_manual(values=c("#F39B7FFF", "#4DBBD5FF", "#00A087FF","#3C5488FF","#E64B35FF","#8491B4FF")) + 
    xlab("Binding score with reference allele") +
    ylab("Binding score with alternative allele") + 
    labs( size = "DeepBind score" ) + 
    xlim(lowerlim,upperlim) + 
    ylim(lowerlim,upperlim) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
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
      complete = T) + 
    guides(colour = guide_legend(override.aes = list(size=7)))
  print(p)
  dev.off()
  
  pdf(paste(fullexonID,"_withouttext_",version,".pdf",sep=""),height=hi,width=wi)
  p=ggplot(data, aes(x=DBscore_WT, y=DBscore_MUT, size=abs(DBscore), color=SNP)) +
    geom_point(alpha=0.7) +
    #scale_color_npg() + 
    scale_color_manual(values=c("#F39B7FFF", "#4DBBD5FF", "#00A087FF","#3C5488FF","#E64B35FF","#8491B4FF")) + 
    xlab("Binding score with reference allele") +
    ylab("Binding score with alternative allele") + 
    labs( size = "DeepBind score" ) +
    xlim(lowerlim,upperlim) + 
    ylim(lowerlim,upperlim) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=12, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
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
      complete = T) + 
    guides(colour = guide_legend(override.aes = list(size=7)))
  print(p)
  dev.off()
}

setwd(outputpath)
bubbleplot(data2plot,outputpath,fullexoninfo,"")
bubbleplot(data2plot0,outputpath,fullexoninfo,"0")
bubbleplot(data2plot00,outputpath,fullexoninfo,"00")








