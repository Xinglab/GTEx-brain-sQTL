#table of content
#1. stacked bar plot of the number of sQTLs in each brain region (SE+A3SS+A5SS)
#2. stacked bar plot of the number of sQTLs significant in 1, 2, 3, ..., 13 brain regions (SE+A3SS+A5SS)
#Of note, the number of sQTL exons here is larger than the number of sQTL exons in the all_info_exon_sQTL_GWAS_disease_in_Brain-XXX_logit_JC_SE_pvalue.txt file of each brain region
#the reason is that the number of sQTL exons here contain sQTL exons with SNPs that cannot be mapped to rsID. These sQTL exons are not included in all_info_exon_sQTL_GWAS_disease_in_Brain-XXX_logit_JC_SE_pvalue.txt but included here
#3. bar plot of the number of unique disease sQTL exons across all regions + heatmap of the number of unique disease sQTL exons in each region + bar plot of the disease enrichment (-log10(pvalue))
#4. bar plot of the number of unique disease sQTL exons across all regions + heatmap of the number of unique disease sQTL exons in each region + bar plot of the disease enrichment (-log10(pvalue))
#5. heat map of beta value for all sQTLs
#6. comparing the region specificity of sQTL exons related to neurodegenerative diseases and other neuropsychiatric diseases

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
#outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,"based_on_high_confidence_GWAS_SNPs",sep="/")

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
    #summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,splicetype,sep="/")
    summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,"based_on_high_confidence_GWAS_SNPs",splicetype,sep="/")
    
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


###################################################################################################################################################################################################
#4. bar plot of the number of unique disease sQTL SNPs across all regions + heatmap of the number of unique disease sQTL SNPs in each region + bar plot of the disease enrichment (-log10(pvalue))#
###################################################################################################################################################################################################
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
#outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,"based_on_high_confidence_GWAS_SNPs",sep="/")

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
    #summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,splicetype,sep="/")
    summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary",PSItype,counttype,"based_on_high_confidence_GWAS_SNPs",splicetype,sep="/")
    
    setwd(summaryinput)
    diseaseexon=read.table(paste("updated_GWAS_loci_in_disease_each_BR_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
    rownames(diseaseexon)=diseaselist
    #remove first two columns and the row of HD
    diseaseexon=diseaseexon[-which(rownames(diseaseexon)=="HD"),]
    uniqueexonnum=diseaseexon[,"unique_disease_GWAS_updated"]
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
  colnames(df.DSEpreregion)=c("Brain.region","Disease","Num.of.sQTL.SNP")
  cellnote=df.DSEpreregion$Num.of.sQTL.SNP    #values of each cell in the heatmap
  df.DSEpreregion$Brain.region <- factor(df.DSEpreregion$Brain.region, levels = rev(formalbrainregionlist))     #fix the order of brain regions
  
  cellnote[cellnote==0]=NA    #we don't plot 0
  hm <- ggplot(data = df.DSEpreregion, aes(x = Disease, y = Brain.region, fill = Num.of.sQTL.SNP)) + geom_tile() + 
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
  colnames(df.DSEtable) <- c("Var1", "Disease","Num.of.sQTL.SNP")
  df.DSEtable$Disease <- factor(df.DSEtable$Disease, levels = df.DSEtable$Disease)   #fix the order of diseases
  barheight.x1=df.DSEtable$Num.of.sQTL.SNP
  barheight.x1[barheight.x1==0]=NA      #if there is no disease related exon, we don't plot it
  bp.x1 <- ggplot(data = df.DSEtable, aes(x = factor(Disease), y = Num.of.sQTL.SNP)) + 
    geom_bar(stat = "identity", aes(fill = Num.of.sQTL.SNP)) + theme_gray() + 
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
  pdf(paste("disease_sQTL_SNP_and_enrichment_",type,".pdf",sep=""),width = 14, height = 14)
  #combine plots
  tmp <- ggplot_gtable(ggplot_build(hm))   #we need to put this part here because this command will call the drawing function. If we don't put it here (don't give a place for output), it will cause error
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  grob.title <- textGrob("", hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 20))
  grid.arrange(bp.x1, hm.clean, bp.x2, nrow = 3, ncol = 1, 
               widths = c(40), heights = c(50, 60, 40), top = grob.title)
  dev.off()
}



#########################################
#5. heat map of beta value for all sQTLs#
#########################################
PSItype="logit"
counttype="JC"
splicetypelist=c("SE","A5SS","A3SS")
typelist=c("pvalue","permutation")

library(gplots)
library(ggplot2)
library(RColorBrewer)
require(scales)

source("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/heatmap.3.R")

vectorsplit=function(x,del,num){
  return(strsplit(x,split=del)[[1]][num])
}

sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/",PSItype,"/",counttype,sep="")
summaryrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/",PSItype,"/",counttype,sep="")
rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot/",PSItype,"/",counttype,sep="")

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
  #get all the sQTL pairs 
  sQTLpair=sQTLdispair=matrix(,nrow=0,ncol=3)
  
  for (splicetype in splicetypelist){
    inputpath=paste(summaryrootinput,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      diseasetemp=subset(temp,temp[,"disease_related"]==1)
      
      temptemp=cbind(temp[,c("exon_short_ID","SNP_info")],rep(splicetype,dim(temp)[1]))
      sQTLpair=rbind(sQTLpair,temptemp)
      
      if (dim(diseasetemp)[1]>0){      #if we have disease related pairs
        diseasetemptemp=cbind(diseasetemp[,c("exon_short_ID","SNP_info")],rep(splicetype,dim(diseasetemp)[1]))
        sQTLdispair=rbind(sQTLdispair,diseasetemptemp)
      }
    }
  }
  sQTLpair=unique(sQTLpair)
  sQTLdispair=unique(sQTLdispair)
  colnames(sQTLpair)=colnames(sQTLdispair)=c("exon","snp","splicetype")
  rownames(sQTLpair)=paste(sQTLpair[,"exon"],sQTLpair[,"snp"],sQTLpair[,"splicetype"],sep="~")
  rownames(sQTLdispair)=paste(sQTLdispair[,"exon"],sQTLdispair[,"snp"],sQTLdispair[,"splicetype"],sep="~")
  
  #get the beta value for all those unique pairs
  betamatrix=matrix(NA,dim(sQTLpair)[1],length(brainregionlist))
  rownames(betamatrix)=rownames(sQTLpair)
  colnames(betamatrix)=brainregionlist
  
  for (i in 1:dim(sQTLpair)[1]){
    if ((i %% 100)==0){
      print(i)
    }
    exonshortID=as.character(sQTLpair[i,"exon"])
    snpID=as.character(sQTLpair[i,"snp"])
    currentST=as.character(sQTLpair[i,"splicetype"])
    #print(i)
    for (br in 1:length(brainregionlist)){
      currentBR=brainregionlist[br]
      #read in the exon result in the current brain region
      inputpath=paste(sqtlrootinput,currentST,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/")
      setwd(inputpath)
      exonresult=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
      row=which(exonresult[,"SNPID"] %in% snpID)[1]
      betamatrix[i,br]=exonresult[row,"Beta"]
    }
  }
  
  disbetamatrix=betamatrix[rownames(sQTLdispair),]
  
  setwd(rootoutput)
  write.table(betamatrix,paste("betamatrix_all_splicetype_",type,".txt",sep=""),sep="\t")
  write.table(disbetamatrix,paste("disease_betamatrix_all_splicetype_",type,".txt",sep=""),sep="\t")
  
  setwd(rootoutput)
  betamatrix=as.matrix(read.table(paste("betamatrix_all_splicetype_",type,".txt",sep=""),sep="\t",header=T))
  disbetamatrix=as.matrix(read.table(paste("disease_betamatrix_all_splicetype_",type,".txt",sep=""),sep="\t",header=T))
  
  #generate heatmap for those pairs
  colnames(betamatrix)=colnames(disbetamatrix)=formalbrainregionlist
  
  setwd(rootoutput)
  #clab=as.character(sapply(rownames(betamatrix),vectorsplit,"~",3))
  #clab[which(clab=="SE")]="red"
  #clab[which(clab=="A3SS")]="blue"
  #clab[which(clab=="A5SS")]="green"
  #clab=t(as.matrix(clab))
  hc <- hclust(dist(betamatrix), method='ave')
  order <- hc$order
  orderedbetamatrix=betamatrix[order,]
  #orderedclab=t(as.matrix(clab[,order]))
  setwd(rootoutput)
  pdf(paste("Region_dependent_original_sQTL_beta_clustering_heatmap_all_splice_type_",type,"_original_beta.pdf",sep=""),height=6,width=6)
  ownbreak = c(seq(range(as.numeric(orderedbetamatrix))[1],-0.8,length=250),seq(-0.79,0.79,length=501),seq(0.8,range(as.numeric(orderedbetamatrix))[2],length=250))
  h=heatmap.3(orderedbetamatrix,                                                 #h here contains the row and column order (h$rowInd and h$colInd) after the clustering. We need this for further plot
              col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(1000), 
              key = TRUE,
              Rowv = FALSE,
              breaks = ownbreak,
              #RowSideColors = orderedclab,
              labRow = NA)
  dev.off()
  
  
  setwd(rootoutput)
  pdf(paste("Region_dependent_original_disease_sQTL_beta_clustering_heatmap_all_splice_type_",type,"_original_beta.pdf",sep=""),height=6,width=6)
  ownbreak = c(seq(range(as.numeric(disbetamatrix))[1],-0.8,length=250),seq(-0.79,0.79,length=501),seq(0.8,range(as.numeric(disbetamatrix))[2],length=250))
  h=heatmap.3(disbetamatrix,                                                 #h here contains the row and column order (h$rowInd and h$colInd) after the clustering. We need this for further plot
              col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(1000), 
              key = TRUE,
              breaks = ownbreak,
              labRow = NA)
  dev.off()
}


#############################################################################################################################
#6. comparing the region specificity of sQTL exons related to neurodegenerative diseases and other neuropsychiatric diseases#
#############################################################################################################################
splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")
counttype="JC"
PSItype="logit"
outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/plot",PSItype,counttype,sep="/")
summaryinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary"

ND=toupper(c("Alzheimer","Amyotrophic lateral sclerosis","Parkinson","frontotemporal dementia","multiple sclerosis"))       #change all the names to upper case for comparison
NPSY=toupper(c("epilepsy","autism","schizophrenia","bipolar","depression",
               "attention deficit hyperactivity disorder","glio","multiple sclerosis","narcolepsy","stroke"))

for (type in typelist){
  NDresult=NPSYresult=matrix(,nrow=0,ncol=4)
  colnames(NDresult)=colnames(NPSYresult)=c("exon","num.sig.region","disease","splicetype")
  for (splicetype in splicetypelist){
    inputpath=paste(summaryinput,PSItype,counttype,splicetype,sep="/")
    setwd(inputpath)
    
    #add the result of neurodegenerative diseases
    for (NDdisease in ND){
      NDdisresult=try(suppressMessages(read.table(paste(gsub("\\s", "_", NDdisease),"_brain_region_specific_exon","_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t")),silent=TRUE)
      if (!(inherits(NDdisresult,"try-error"))) {       #if we have this result
        temp=cbind(rownames(NDdisresult),apply(NDdisresult,1,sum),rep(NDdisease,dim(NDdisresult)[1]),rep(splicetype,dim(NDdisresult)[1]))
        NDresult=rbind(NDresult,temp)
      }
    }
    
    #add the result of other neuropsychiatric diseases
    for (NPSYdisease in NPSY){
      NPSYdisresult=try(suppressMessages(read.table(paste(gsub("\\s", "_", NPSYdisease),"_brain_region_specific_exon","_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t")),silent=TRUE)
      if (!(inherits(NPSYdisresult,"try-error"))) {       #if we have this result
        temp=cbind(rownames(NPSYdisresult),apply(NPSYdisresult,1,sum),rep(NPSYdisease,dim(NPSYdisresult)[1]),rep(splicetype,dim(NPSYdisresult)[1]))
        NPSYresult=rbind(NPSYresult,temp)
      }
    }
  }
  nd=sort(as.numeric(NDresult[,"num.sig.region"]))
  npsy=sort(as.numeric(NPSYresult[,"num.sig.region"]))
  p.value=t.test(nd,npsy)$p.value
  print(paste(type,p.value))
}
#[1] "pvalue 0.618722252994062"
#[1] "permutation 0.517855848952797"
#Conclusion: the region specificity of neurodegenerative diseases related sQTL exons is not significantly different from other neuropsychiatric diseases


#disease classification:
#ND: AD, frontotemporal dementia, Parkinson’s disease (PD), Huntington’s disease (HD), amyotrophic lateral sclerosis (ALS), multiple sclerosis (MS). 
#Other: Epilepsy, Stroke, Narcolepsy, Attention deficit hyperactivity disorder (ADHD), Autism, Schizophrenia, Bipolar, Depression, glioblastoma)
#MS is neurodegenerative:
#(1) Multiple sclerosis is primarily a neurodegenerative disease.
#(2) The role of brain vasculature in neurodegenerative disorders 
#(3) Multiple Sclerosis: An Immune or Neurodegenerative Disorder?
#(4) Unravelling neurodegeneration in multiple sclerosis (The lancet neurology)






