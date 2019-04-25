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

#sQTLexoninputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.1_sQTL_disease_exon_annotation/1_disease_exon_summary/result"
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/sQTL_significance_vs_brain_region_specificity"
summaryinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary"
sQTLinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"

library(ggplot2)

for (splicetype in splicetypelist){
  for (type in typelist){
    #1. get all sQTL exons
    sQTLexon=c()
    inputpath=paste(summaryinput,PSItype,counttype,splicetype,sep="/")
    setwd(inputpath)
    for (i in 1:length(brainregionlist)){
      temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
      sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
    }
    sQTLexon=unique(sQTLexon)
    
    #2. calculate average sQTL significance across all significant regions for each sQTL exon & count the number of significant regions#
    num.sig.region=rep(0,length(sQTLexon))
    exoninfo=matrix(0,length(sQTLexon),length(brainregionlist))     #store the p values
    colnames(exoninfo)=brainregionlist
    rownames(exoninfo)=sQTLexon
    sigexoninfo=exoninfo          #a indicator matrix with 1 means significant and 0 means insignificant
    
    for (e in 1:length(sQTLexon)){
      #print(e)
      for (br in 1:length(brainregionlist)){
        #get the p value cutoff
        if (type=="pvalue"){
          cutoff=10^-5
        }
        if (type=="permutation"){
          setwd(paste(sQTLinput,"/",PSItype,"/",counttype,"/",splicetype,"/",brainregionlist[br],sep=""))
          cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
        }
        
        #print(br)
        shortexonID=paste("SE_",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="")
        path=paste(sQTLinput,PSItype,counttype,splicetype,brainregionlist[br],paste("Glimmps_each_exon_cis_",brainregionlist[br],sep=""),sep="/")
        setwd(path)
        temp=read.table(paste(shortexonID,".asso",sep=""),header=T)
        p.value=min(temp[,"pvals.lm"])
        exoninfo[e,br]=p.value
        if (p.value<=cutoff){
          num.sig.region[e]=num.sig.region[e]+1
          sigexoninfo[e,br]=1
        }
      }
    }
    
    significance_mean=rep(NA,length(sQTLexon))
    significance_median=rep(NA,length(sQTLexon))
    significance_min=rep(NA,length(sQTLexon))
    significance_max=rep(NA,length(sQTLexon))
    for (i in 1:length(sQTLexon)){
      pvaluelist=exoninfo[i,which(sigexoninfo[i,]==1)]
      significance_mean[i]=mean(pvaluelist)
      significance_median[i]=median(pvaluelist)
      significance_min[i]=min(pvaluelist)
      significance_max[i]=max(pvaluelist)
    }
    significance_mean=-log10(significance_mean)
    significance_median=-log10(significance_median)
    significance_min=-log10(significance_min)
    significance_max=-log10(significance_max)
    
    #4. make the plot
    for (label in c("mean","median","min","max")){
      if (label=="mean"){
        significance=significance_mean
      }
      if (label=="median"){
        significance=significance_median
      }
      if (label=="min"){
        significance=significance_min
      }
      if (label=="max"){
        significance=significance_max
      }
      data2plot=cbind(significance,num.sig.region)
      df=as.data.frame(data2plot)
      df$significance=as.numeric(as.matrix(df$significance))
      df$num.sig.region <- factor(df$num.sig.region, levels = as.character(seq(1:13)))   #fix the order of diseases
      
      setwd(outputpath)
      
      ###bar plot###
      pdf(paste(splicetype,"_",type,"_",label,"_","sQTL_significance_vs_brain_region_specificity_bar_plot.pdf",sep=""),height=6,width=6)
      p <- ggplot(df, aes(x=num.sig.region, y=significance, fill=num.sig.region)) + 
        geom_boxplot() + 
        labs(y = expression('sQTL significance (-log'[10]*'(P))'), x="number of significant regions") + 
        ggtitle(paste(splicetype,type,label)) +
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
      
      ###violin plot###
      pdf(paste(splicetype,"_",type,"_",label,"_","sQTL_significance_vs_brain_region_specificity_violin_plot.pdf",sep=""),height=6,width=6)
      p <- ggplot(df, aes(x=num.sig.region, y=significance, fill=num.sig.region)) + 
        geom_violin() + stat_summary(fun.y=mean, geom="point", shape=23, size=2) + geom_boxplot(width=0.1) +
        labs(y = expression('sQTL significance (-log'[10]*'(P))'), x="number of significant regions") + 
        ggtitle(paste(splicetype,type,label)) +
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
}

