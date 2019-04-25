outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/10.4_radar_plot_for_RBP_expression_and_sQTL_significance/result"

library(NMF)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(fmsb)

splicetype="SE"
sqtlinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                splicetype,sep="/")

#exon ID conversion table
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/SE/pvalue")
exonIDconversion=read.table("SE_pvalue_short_long_exon_ID_conversion.txt",sep="\t")
rownames(exonIDconversion)=exonIDconversion[,1]

RBP=c("ENSG00000100320.18_RBFOX2","ENSG00000139910.15_NOVA1","ENSG00000165119.14_HNRNPK")
exonsnp=c("SE_189436~17_44051924_G_A_b37",
          "SE_189436~17_44051612_A_G_b37",
          "SE_303254~12_51593616_T_G_b37",
          "SE_108190~2_242793_G_C_b37",    #SH3YL1 on splice site
          "SE_108171~2_224919_A_G_b37",    #SH3YL1 on splice site
          "SE_229974~3_113012797_G_A_b37",  #CFAP44 on dinucleotide
          "SE_222152~6_30708695_C_T_b37",   #flot1 & one significant SNP within 300bp
          "SE_344917~11_47434986_G_A_b37", #SLC39A13 & one significant SNP within 300bp
          "SE_128125~5_139059562_A_G_b37",  #CXXC5
          "SE_20599~3_39540727_C_A_b37",    #MOBP
          "SE_263001~6_30172513_A_C_b37",   #TRIM26
          "SE_359526~20_62320968_T_C_b37",  #RTEL1
          "SE_9912~5_176798306_G_A_b37",    #RGS14
          "SE_222152~6_30644137_T_C_b37",   #flot1 & sQTL SNP
          "SE_344917~11_47509137_G_C_b37")  #SLC39A13 & sQTL SNP


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

####################
#read in annotation#
####################
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


#################################
#Spider plot for gene expression#
#################################
for (r in 1:length(RBP)){
  currentRBP=RBP[r]
  
  data2plot=as.data.frame(matrix(0 , nrow=length(currentRBP) , ncol=length(formalbrainregionlist)))
  rownames(data2plot)=currentRBP
  colnames(data2plot)=formalbrainregionlist
  
  for (b in 1:length(formalbrainregionlist)){
    currentBR=formalbrainregionlist[b]
    exprRBP=edata[currentRBP,which(brainregion %in% currentBR)]
    data2plot[currentRBP,b]=mean(exprRBP,na.rm=T)
  }
  
  data2plot=rbind(rep(ceiling(max(data2plot)),length(formalbrainregionlist)) , rep(0,length(formalbrainregionlist)) , data2plot)
  
  setwd(outputpath)
  write.table(data2plot,paste("RBPexpr~",currentRBP,".txt",sep=""),sep="\t",row.names=F)
  
  setwd(outputpath)
  pdf(paste("RBPexpr~",currentRBP,".pdf",sep=""),height=6,width=6)
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.7,0.5,0.1,0.9) )      #we only use one color for each plot so the second one is not used
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.7,0.5,0.1,0.4) )
  radarchart( data2plot  , axistype=1 , 
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(from=0,to=ceiling(max(data2plot)),length.out=5), cglwd=4,
              #custom labels
              vlcex=1,
              centerzero=T
  )
  #legend(x=1, y=1, legend = rownames(data2plot), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
  dev.off()
}


###################################
#Spider plot for sQTL significance#
###################################
for (es in 1:length(exonsnp)){
  currentpair=exonsnp[es]
  currentexon=strsplit(currentpair,split="~")[[1]][1]
  currentexonsymbol=strsplit(as.character(exonIDconversion[currentexon,2]),split="\\|")[[1]][3]
  currentsnp=strsplit(currentpair,split="~")[[1]][2]
  
  data2plot=as.data.frame(matrix(0 , nrow=length(currentpair) , ncol=length(formalbrainregionlist)))
  rownames(data2plot)=exonsnp[es]
  colnames(data2plot)=formalbrainregionlist
  
  for (b in 1:length(formalbrainregionlist)){
    currentBR=brainregionlist[b]
    setwd(paste(sqtlinput,currentBR,paste("Glimmps_each_exon_cis_",currentBR,sep=""),sep="/"))
    exonresult=read.table(paste(currentexon,".asso",sep=""),sep="\t",header=T)
    row=which(exonresult[,"SNPID"] %in% currentsnp)[1]
    p.value=exonresult[row,"pvals.lm"]
    
    data2plot[currentpair,b]=p.value
  }
  
  data2plot=-log10(data2plot)
  data2plot=rbind(rep(ceiling(max(data2plot)),length(formalbrainregionlist)) , rep(0,length(formalbrainregionlist)) , data2plot)
  
  setwd(outputpath)
  write.table(data2plot,paste("sQTLsig~",currentpair,"_",currentexonsymbol,".txt",sep=""),sep="\t",row.names=F)
  
  numseg=ceiling(max(data2plot))     #number of segments for each axis (default 4). We want to show the significance cutoff so we need to change this for different events
  
  setwd(outputpath)
  pdf(paste("sQTLsig~",currentpair,"_",currentexonsymbol,".pdf",sep=""),height=6,width=6)
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.7,0.5,0.1,0.9) )
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.7,0.5,0.1,0.4) )
  radarchart( data2plot  , axistype=1 , 
              maxmin = T,
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="black", 
              caxislabels=seq(from=0,to=ceiling(max(data2plot)),length.out=numseg+1),        
              cglwd=4,
              #custom labels
              vlcex=1,
              #other stuff
              centerzero=T,
              seg=numseg,
              title=paste(currentexonsymbol,currentexon,currentsnp,sep="~")
  )
  dev.off()
}

