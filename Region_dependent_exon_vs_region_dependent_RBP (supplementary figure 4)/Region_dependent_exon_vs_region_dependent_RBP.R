#The purpose of this code is to plot the gene expression/PSI value of genes/exons across different brain regions.
#Here we use uncorrected PSI and gene expression because the result of this analysis is used before any correction.
#The result here is just a sanity check of our data. 
outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/10.1_region_dependent_exon_vs_region_dependent_RBP/result"
library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)

#######################################################################################################
#########################
#read in gene expression#
#########################
exprinputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM"
expressiontype="allSF"
setwd(exprinputpath)
if (expressiontype=="allgene"){
  expr=read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",header=T,sep="\t",check.names=F)
}
if (expressiontype=="allRBP"){
  expr=read.table("all_RBP_normalized_expression.txt",header=T,sep="\t",check.names=F)
}
if (expressiontype=="allSF"){
  expr=read.table("SF_normalized_expression.txt",header=T,sep="\t",check.names=F)
}
rownames(expr)=paste(expr[,"Name"],expr[,"Description"],sep="_")
edata = as.matrix(expr[,-c(1,2)])
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]
columnname=rep(NA,dim(edata)[2])
for (i in 1:dim(edata)[2]){
  columnname[i]=as.character(IDconversion[which(IDconversion[,"sampleID"] %in% colnames(edata)[i]),"SRRlist"])
}
colnames(edata)=columnname

###################
#read in phenotype#
###################
phenotypepath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"
setwd(phenotypepath)
orisampleID_SRRID=read.table("sampleID_SRRID_brain.txt",sep="\t",header=T)
oribrainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)

#get the phenotype information of samples in the expr matrix (not all 1409 sampes are in the expr data)
brainregion=oribrainregion[,as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"])]
sampleID_SRRID=orisampleID_SRRID[as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"]),]

#change brain region labels
brainregion[which(brainregion=="Brain - Cerebellar Hemisphere")]="Cerebellar Hemisphere"
brainregion[which(brainregion=="Brain - Cerebellum")]="Cerebellum"
brainregion[which(brainregion=="Brain - Cortex")]="Cortex"
brainregion[which(brainregion=="Brain - Frontal Cortex (BA9)")]="Frontal Cortex BA9"
brainregion[which(brainregion=="Brain - Anterior cingulate cortex (BA24)")]="Anterior cingulate cortex BA24"
brainregion[which(brainregion=="Brain - Hypothalamus")]="Hypothalamus"
brainregion[which(brainregion=="Brain - Caudate (basal ganglia)")]="Caudate basal ganglia"
brainregion[which(brainregion=="Brain - Amygdala")]="Amygdala"
brainregion[which(brainregion=="Brain - Spinal cord (cervical c-1)")]="Spinal cord cervical c-1"
brainregion[which(brainregion=="Brain - Nucleus accumbens (basal ganglia)")]="Nucleus accumbens basal ganglia"
brainregion[which(brainregion=="Brain - Substantia nigra")]="Substantia nigra"
brainregion[which(brainregion=="Brain - Putamen (basal ganglia)")]="Putamen basal ganglia"
brainregion[which(brainregion=="Brain - Hippocampus")]="Hippocampus"

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
#############
#read in PSI#
#############
splicetype="SE"
counttype="JC"
PSIinputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.2_data_inspection_imputation/result/",counttype,sep="")
setwd(PSIinputpath)
PSI=read.table(paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t",header=T)
PSI=PSI[,colnames(edata)]

######
#plot#
######
cbPalette <- c("#E53537",     #GTEx Artery color/Amygdala
               "#EEA760",     #GTEx adipose tissue color/ACC
               "#C9E5C3",     #GTEx pituitary color/Caudate
               "#9ACB3C",     #GTEx lung color/CH
               "#028B45",     #GTEx thyroid color/Cerebellum
               "#FFD923",     #GTEx nerve color/cortex
               "#F0EC68",     #GTEx brain color/frontal cortex
               "#6D67AF",     #GTEx muscle color/hippocampus
               "#7F388D",     #GTEx Heart color/Hypothalamus
               "#2999D5",     #GTEx Skin color/Nucleus accumbens
               "#4EC3C7",     #GTEx breast color/Putamen
               "#CD8ABC",     #GTEx LCL color/spinal cord
               "#F287B7")     #GTEx blood color/Substantia nigra

gene2plot=c("ENSG00000121774.13_KHDRBS1",      #Sam68
            "ENSG00000131773.9_KHDRBS3")       #T-STAR       

exon2plot=c("135968|ENSG00000021645.18_3|NRXN3|chr14|+|80158515|80158605|80130120|80130292|80163972|80164280",     #this is the NRXN3 in the paper
            "30935|ENSG00000110076.18_2|NRXN2|chr11|-|64393934|64394024|64390224|64390550|64397873|64398045",      #this is the NRXN2 in the paper
            "30981|ENSG00000110076.18_2|NRXN2|chr11|-|64421167|64421194|64419506|64419626|64427803|64428007",
            "30988|ENSG00000110076.18_2|NRXN2|chr11|-|64444464|64444509|64435914|64436076|64453117|64453419",
            "30991|ENSG00000110076.18_2|NRXN2|chr11|-|64444485|64444509|64435914|64436076|64453117|64453419",
            "31005|ENSG00000110076.18_2|NRXN2|chr11|-|64457876|64457948|64453117|64453419|64460318|64460348",
            "31006|ENSG00000110076.18_2|NRXN2|chr11|-|64457876|64457948|64453117|64453419|64465246|64465264",
            "31010|ENSG00000110076.18_2|NRXN2|chr11|-|64460318|64460348|64457876|64457948|64465246|64465264",
            "242960|ENSG00000179915.22_3|NRXN1|chr2|-|50282092|50282182|50280408|50280728|50318460|50318632",     #this is the NRXN1 in the paper
            "242970|ENSG00000179915.22_3|NRXN1|chr2|-|50755762|50755789|50733632|50733755|50758364|50758568")

for (gene in gene2plot){
  geneexpr=edata["ENSG00000131773.9_KHDRBS3",]/edata["ENSG00000121774.13_KHDRBS1",]      #here we plot the ration instead of expression level
  #geneexpr=edata[gene,]
  data2plot=t(rbind(brainregion,geneexpr))
  colnames(data2plot)=c("Brain.region","Expression")
  colnames(data2plot)=c("Brain.region","Expression")
  data2plot=as.data.frame(data2plot)
  df=melt(data2plot)
  df$Brain.region=factor(df$Brain.region,levels=formalbrainregionlist)
  df$Expression=as.numeric(as.character(df$Expression))
  
  #plot in pdf#
  setwd(outputpath)
  #pdf(paste("Region_dependent_gene_expression_of_",gene,".pdf",sep=""),height=8,width=8)
  pdf(paste("Region_dependent_gene_expression_of_","ratio",".pdf",sep=""),height=8,width=8)
  limit=max(abs(floor(range(as.numeric(as.character(df[,"Expression"]))))[1]),abs(ceiling(range(as.numeric(as.character(df[,"Expression"]))))[2]))
  p <- ggplot(df, aes(x=Brain.region, y=Expression, fill=Brain.region, color=Brain.region)) + 
    geom_boxplot() + 
    scale_fill_manual(values=cbPalette) + 
    scale_color_manual(values=cbPalette) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) + geom_boxplot(width=0.2) + geom_jitter(shape=16, position=position_jitter(0.05)) + 
    #labs(x = "Brain region", y="Expression") + 
    labs(x = "Brain region", y="Expression Ratio") + 
    scale_y_continuous(breaks = round(seq(-limit,limit, by = 0.5),1)) + 
    ylim(0,2.5) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=14, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 1, vjust = 1, size=14, angle=60, face="bold"),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14, face="bold"),
      axis.line = element_line(colour = "black"),
      #background
      panel.background = element_blank(),
      #legend
      #legend.position="none",
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)+
    #ggtitle(gene)
    ggtitle("Ratio of T-STAR:Sam68")
  print(p)
  dev.off()
}


for (exon in exon2plot){
  data2plot=t(rbind(brainregion,PSI[exon,]))
  colnames(data2plot)=c("Brain.region","PSI")
  data2plot=as.data.frame(data2plot)
  df=melt(data2plot)
  df$Brain.region=factor(df$Brain.region,levels=formalbrainregionlist)
  df$PSI=as.numeric(as.character(df$PSI))
  
  #plot in pdf#
  setwd(outputpath)
  pdf(paste("Region_dependent_alternative_splicing_of_",exon,".pdf",sep=""),height=8,width=8)
  limit=max(abs(round(range(as.numeric(as.character(df[,"PSI"]))))[1]),abs(round(range(as.numeric(as.character(df[,"PSI"]))))[2]))
  p <- ggplot(df, aes(x=Brain.region, y=PSI, fill=Brain.region, color=Brain.region)) + 
    geom_boxplot() + 
    scale_fill_manual(values=cbPalette) + 
    scale_color_manual(values=cbPalette) + 
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) + geom_boxplot(width=0.2) + geom_jitter(shape=16, position=position_jitter(0.05)) + 
    labs(x = "Brain region", y="PSI") + 
    #scale_y_continuous(breaks = round(seq(-limit,limit, by = 4),1)) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=14, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=14, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 1, vjust = 1, size=14, angle=60, face="bold"),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14, face="bold"),
      axis.line = element_line(colour = "black"),
      #background
      panel.background = element_blank(),
      #legend
      #legend.position="none",
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)+
    ggtitle(exon)
  print(p)
  dev.off()
}


