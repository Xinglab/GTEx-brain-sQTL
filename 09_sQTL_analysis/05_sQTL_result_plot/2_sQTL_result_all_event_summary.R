#table of content
#1. stacked bar plot of the number of sQTLs in each brain region (SE+A3SS+A5SS)
#2. stacked bar plot of the number of sQTLs significant in 1, 2, 3, ..., 13 brain regions (SE+A3SS+A5SS)
#3. bar plot of the number of unique disease sQTL exons across all regions + heatmap of the number of unique disease sQTL exons in each region + bar plot of the disease enrichment (-log10(pvalue))

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
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")

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
    summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,splicetype,sep="/")
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
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")
rootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run",PSItype,counttype,sep="/")
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

#####################################################################################################################################################################################################
#3. bar plot of the number of unique disease sQTL exons across all regions + heatmap of the number of unique disease sQTL exons in each region + bar plot of the disease enrichment (-log10(pvalue))#
#####################################################################################################################################################################################################
library("data.table")
library(gaston)
library(boot)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
#library(cowplot)

splicetypelist=c("SE","A5SS","A3SS")
typelist=c("pvalue","permutation")
counttype="JC"
PSItype="logit"
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")

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

diseaselist=c("AD","ALS","PD","FTD","HD","Epilepsy","Autism","Schizophrenia","Bipolar","Depression",
              "ADHD","Glioma/Glioblastoma","MS","Narcolepsy","Stroke")
columnname=c(paste(splicetypelist[1],diseaselist[-which(diseaselist=="HD")],sep="_"),
             paste(splicetypelist[2],diseaselist[-which(diseaselist=="HD")],sep="_"),
             paste(splicetypelist[3],diseaselist[-which(diseaselist=="HD")],sep="_"))


for (type in typelist){
  #1. the number of unique disease sQTL exons across all regions in SE, A3SS, and A5SS
  DSEtable=matrix(0,length(splicetypelist),length(diseaselist)-1)
  rownames(DSEtable)=splicetypelist
  colnames(DSEtable)=diseaselist[-which(diseaselist=="HD")]
  #2. the number of unique disease sQTL exons in each region
  DSEperregion=matrix(NA,length(brainregionlist),0)
  
  for (splicetype in splicetypelist){
    summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,splicetype,sep="/")
    setwd(summaryinput)
    diseaseexon=read.table(paste("upadted_exon_in_disease_each_BR_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
    rownames(diseaseexon)=diseaselist
    #remove first two columns and the row of HD
    diseaseexon=diseaseexon[-which(rownames(diseaseexon)=="HD"),]
    uniqueexonnum=diseaseexon[,"unique_disease_exon_updated"]
    diseaseexon=diseaseexon[,-c(1,2)]
    colnames(diseaseexon)=formalbrainregionlist
    
    DSEtable[splicetype,]=uniqueexonnum
    DSEperregion=cbind(DSEperregion,t(diseaseexon))
  }
  DSEtable=cbind(t(DSEtable["SE",]),t(DSEtable["A5SS",]),t(DSEtable["A3SS",]))
  
  #3. disease enrichment p value for all diseases
  setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/disease_GWAS_enrichment")
  diseaseenrichment=read.table(paste("disease_enrichment_pvalue_",type,".txt",sep=""),sep="\t")
  colnames(diseaseenrichment)=diseaselist[-which(diseaselist=="HD")]
  diseaseenrichment=cbind(diseaseenrichment["SE",],diseaseenrichment["A5SS",],diseaseenrichment["A3SS",])
  rownames(diseaseenrichment)=NULL
  
  colnames(DSEtable)=columnname
  colnames(DSEperregion)=columnname
  colnames(diseaseenrichment)=columnname
  
  #make plot (the order is SE, A3SS, A5SS)
  #prepare heatmap
  df.DSEpreregion=melt(DSEperregion)
  colnames(df.DSEpreregion)=c("Brain.region","Disease","Num.of.sQTL.exon")
  cellnote=df.DSEpreregion$Num.of.sQTL.exon    #values of each cell in the heatmap
  df.DSEpreregion$Brain.region <- factor(df.DSEpreregion$Brain.region, levels = rev(formalbrainregionlist))     #fix the order of brain regions
  
  cellnote[cellnote==0]=NA    #we don't plot 0
  hm <- ggplot(data = df.DSEpreregion, aes(x = Disease, y = Brain.region, fill = Num.of.sQTL.exon)) + geom_tile() + 
    scale_fill_distiller(name = "",      #no legend name
                         palette = "Reds", direction = 1, na.value = "transparent") +
    scale_x_discrete(breaks = unique(df.DSEpreregion$Disease), labels = unique(df.DSEpreregion$Disease)) + theme_gray() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 15), 
          legend.key.size = unit(1,"cm"),
          #legend.text = element_text(size = 7), 
          legend.text = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1),
          #axis.text.y = element_text(angle = 30, hjust = 1),
          text = element_text(size=18, face="bold"), 
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5)) + 
    geom_text(aes(label = cellnote, fontface="bold"))     #add number of exons onto each cell
  # Remove legend and background from heatmap
  hm.clean <- hm +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position="none")
  
  #Prepare x axis barplot 1: disease sQTL exon across all regions
  df.DSEtable=melt(DSEtable)
    #data.frame(rbind(colnames(DSEtable),DSEtable))
  colnames(df.DSEtable) <- c("Var1", "Disease","Num.of.sQTL.exon")
  df.DSEtable$Disease <- factor(df.DSEtable$Disease, levels = df.DSEtable$Disease)   #fix the order of diseases
  barheight.x1=df.DSEtable$Num.of.sQTL.exon
  barheight.x1[barheight.x1==0]=NA      #if there is no disease related exon, we don't plot it
  bp.x1 <- ggplot(data = df.DSEtable, aes(x = factor(Disease), y = Num.of.sQTL.exon)) + 
    geom_bar(stat = "identity", aes(fill = Num.of.sQTL.exon)) + theme_gray() + 
    geom_text(aes(label=barheight.x1, fontface="bold", size=3), vjust=1.5) +    #add bar height onto each bar
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.text.x = element_text(size = 14, face="bold", angle=60, hjust=0),     #name of each disease
          #axis.text.x = element_blank(),
          #axis.title.x = element_text(size = 14, margin = margin(10,0,0,0)),
          axis.title.x = element_blank(),
          legend.position = "none") +
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) + 
    scale_x_discrete(position="top") +   #put the axis on top instead of bottom
    labs(x = "Disease") 
  #bp.x.flip <- switch_axis_position(bp.x, "x")   #flip x axis
  
  #Prepare x axis barplot 2: disease enrichment
  df.diseaseenrichment=melt(diseaseenrichment)
  #data.frame(rbind(colnames(DSEtable),DSEtable))
  colnames(df.diseaseenrichment) <- c("Disease","transformed.p.value")
  df.diseaseenrichment$Disease <- factor(df.diseaseenrichment$Disease, levels = df.diseaseenrichment$Disease)   #fix the order of diseases
  df.diseaseenrichment$transformed.p.value <- -log10(df.diseaseenrichment$transformed.p.value)
  barheight.x2=df.diseaseenrichment$transformed.p.value
  barheight.x2[barheight.x2==0]=NA      #if there is no disease related exon, we don't plot it
  bp.x2 <- ggplot(data = df.diseaseenrichment, aes(x = factor(Disease), y = transformed.p.value)) + 
    geom_bar(stat = "identity", aes(fill = transformed.p.value)) + theme_gray() + 
    #geom_text(aes(label=barheight.x2, fontface="bold", size=3), vjust=1.5) +    #add bar height onto each bar
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          #axis.text.x = element_text(size = 16, face="bold", angle=60, hjust=0),     #name of each disease
          axis.text.x = element_blank(),
          #axis.title.x = element_text(size = 14, margin = margin(10,0,0,0)),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) + 
    scale_x_discrete(position="top") +   #put the axis on top instead of bottom
    labs(x = "Disease") + 
    geom_hline(yintercept = -log10(0.05)) +  #add significance cutoff line
    scale_y_reverse()
  #bp.x.flip <- switch_axis_position(bp.x, "x")   #flip x axis
  
  setwd(outputpath)
  pdf(paste("disease_sQTL_exon_and_enrichment_",type,".pdf",sep=""),width = 14, height = 14)
  #combine plots
  tmp <- ggplot_gtable(ggplot_build(hm))   #we need to put this part here because this command will call the drawing function. If we don't put it here (don't give a place for output), it will cause error
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  grob.title <- textGrob("", hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
  grid.arrange(bp.x1, hm.clean, bp.x2, nrow = 3, ncol = 1, 
               widths = c(40), heights = c(50, 60, 40), top = grob.title)
  dev.off()
}






