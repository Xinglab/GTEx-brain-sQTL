#The purpose of this code:
#we group the samples into two groups: samples from regions where the sQTL is significant and samples from regions where the sQTL is insignificant
#then we get the expression of those samples in the two groups and perform a differential expression analysis

splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
type="pvalue"       #pvalue or permutation
PSItype="logit"
counttype="JC"

library(robust)
library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)
library(corrplot)
library(ggsci)

library("data.table")
library(gaston)
library(boot)

'%!in%' <- function(x,y)!('%in%'(x,y))
source("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/heatmap.3.R")

vectorsplit=function(x,del,pos){
  return(strsplit(x,split=del)[[1]][pos])
}

outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_region_dependent_sQTL_RBP_plot_differential_expression_version",
                 splicetype,type,sep="/")
command=paste("mkdir -p",outputpath)
system(command)

sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                    splicetype,sep="/")

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

#######################
#splicing and genotype#
#######################
PSItype="logit"
counttype="JC"
inputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/5_correction_for_technical_confounders/result/",counttype,sep="")
summaryinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/",PSItype,"/",counttype,"/",splicetype,sep="")
rootsqtl="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run"
totalcountinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")
GWASdbpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases"
GWASdbname="gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
genoplinkprefix="Genotype_V7_plink_maf0.05"
LDpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/LD"
sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                    splicetype,sep="/")

#read in the PSI value
inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")
setwd(inputpath)
PSI=read.table(inputPSI,sep="\t",header=T)              #PSI as the normalized inclusion count
setwd(totalcountinput)
totalRC=read.table(inputtotalRC,sep="\t",header=T)
#transform logit PSI back to PSI
if (PSItype=="logit"){
  temp=inv.logit(as.numeric(as.matrix(PSI)))
  originalPSI=matrix(temp,dim(PSI)[1],dim(PSI)[2])
  rownames(originalPSI)=rownames(PSI)
  colnames(originalPSI)=colnames(PSI)
  PSI=as.data.frame(originalPSI)
}

#################
#gene expression#
#################
#read in the joblist
joblistinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(joblistinput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")
#read in the expression table
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM")
normTPM=read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",sep="\t",header=T,check.names=F)
normTPM=as.matrix(normTPM)

#change RBP name (the name of some RBPs in the DeepBind table is different with their name in the gene expression table. They are the same RBP but use different names)
normTPM[which(normTPM[,2] %in% "CELF4"),2]="BRUNOL4"
normTPM[which(normTPM[,2] %in% "ELAVL1"),2]="HuR"
normTPM[which(normTPM[,2] %in% "TUT1"),2]="STAR-PAP"
normTPM[which(normTPM[,2] %in% "CELF5"),2]="BRUNOL5"
normTPM[which(normTPM[,2] %in% "CELF6"),2]="BRUNOL6"
normTPM[which(normTPM[,2] %in% "HNRNPLL"),2]="HNRPLL"

deepbindRBPlist=as.character(unique(uniquejoblist[,"RBP"]))
#gene ID conversion table
RBPIDconversion=matrix(NA,length(deepbindRBPlist),4)
colnames(RBPIDconversion)=c("DeepBindID","GeneSymbol","EnsemblID","edataID")
RBPIDconversion[,"DeepBindID"]=deepbindRBPlist
RBPIDconversion[,"GeneSymbol"]=sapply(deepbindRBPlist,vectorsplit,del="_",2)
Ensembl2symbol=cbind(as.character(normTPM[,1]),as.character(normTPM[,2]),paste(as.character(normTPM[,1]),as.character(normTPM[,2]),sep="_"))
rownames(Ensembl2symbol)=Ensembl2symbol[,2]
RBPIDconversion[,"EnsemblID"]=Ensembl2symbol[RBPIDconversion[,"GeneSymbol"],1]
RBPIDconversion[,"edataID"]=Ensembl2symbol[RBPIDconversion[,"GeneSymbol"],3]
#######################################################################################################################################################################
#add RBFOX2 and RBFOX3 into the table. Although we don't have Deepbind model for RBFOX2/3, they share very similar motif with RBFOX1. So we still want to include them#
#######################################################################################################################################################################
RBFOX2=c("D00210.001_RBFOX1","RBFOX2","ENSG00000100320.18","ENSG00000100320.18_RBFOX2")
RBFOX3=c("D00210.001_RBFOX1","RBFOX3","ENSG00000167281.14","ENSG00000167281.14_RBFOX3")
RBPIDconversion=rbind(rbind(RBPIDconversion,RBFOX2),RBFOX3)
################################################################################
edata=normTPM
rownames(edata)=paste(as.character(normTPM[,1]),as.character(normTPM[,2]),sep="_")
edata=edata[,-c(1,2)]
temp=matrix(as.numeric(edata),dim(edata)[1],dim(edata)[2])
rownames(temp)=rownames(edata)
colnames(temp)=colnames(edata)
edata=temp

##########################
#brain region information#
##########################
phenotypepath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"
setwd(phenotypepath)
oribrainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)
oribrainregion[which(oribrainregion=="Brain - Amygdala")]="Brain-Amygdala"
oribrainregion[which(oribrainregion=="Brain - Anterior cingulate cortex (BA24)")]="Brain-AnteriorcingulatecortexBA24"
oribrainregion[which(oribrainregion=="Brain - Caudate (basal ganglia)")]="Brain-Caudatebasalganglia"
oribrainregion[which(oribrainregion=="Brain - Cerebellar Hemisphere")]="Brain-CerebellarHemisphere"
oribrainregion[which(oribrainregion=="Brain - Cerebellum")]="Brain-Cerebellum"
oribrainregion[which(oribrainregion=="Brain - Cortex")]="Brain-Cortex"
oribrainregion[which(oribrainregion=="Brain - Frontal Cortex (BA9)")]="Brain-FrontalCortexBA9"
oribrainregion[which(oribrainregion=="Brain - Hippocampus")]="Brain-Hippocampus"
oribrainregion[which(oribrainregion=="Brain - Hypothalamus")]="Brain-Hypothalamus"
oribrainregion[which(oribrainregion=="Brain - Nucleus accumbens (basal ganglia)")]="Brain-Nucleusaccumbensbasalganglia"
oribrainregion[which(oribrainregion=="Brain - Putamen (basal ganglia)")]="Brain-Putamenbasalganglia"
oribrainregion[which(oribrainregion=="Brain - Spinal cord (cervical c-1)")]="Brain-Spinalcordcervicalc-1"
oribrainregion[which(oribrainregion=="Brain - Substantia nigra")]="Brain-Substantianigra"
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM")
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

colnames(edata)=as.character(IDconversion[colnames(edata),"SRRlist"])
edataRBP=as.matrix(edata[RBPIDconversion[,"edataID"],])
rownames(edataRBP)=RBPIDconversion[,"edataID"]
brainregion=oribrainregion[,as.character(colnames(edata))]


#get the full exon information
#get all sQTL exons
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/summary/logit/JC"
sQTLexon=c()
inputpath=paste(rootinput,splicetype,sep="/")
setwd(inputpath)
for (i in 1:length(brainregionlist)){
  temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
  sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
}
sQTLexon=unique(sQTLexon)
shortIDlist=rep(NA,length(sQTLexon))
exonsymbol=rep(NA,length(sQTLexon))
for (e in 1:length(sQTLexon)){
  shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
  exonsymbol[e]=strsplit(sQTLexon[e],split="\\|")[[1]][3]
}
exonIDconversion=cbind(sQTLexon,shortIDlist,exonsymbol)
rownames(exonIDconversion)=shortIDlist
colnames(exonIDconversion)=c("fullID","shortID","gene.symbol")

setwd(outputpath)
write.table(RBPIDconversion,"RBP_ID_conversion.txt",sep="\t")
write.table(exonIDconversion,"Exon_ID_conversion.txt",sep="\t")
###################################################################################################################################################
###########
#section 1#
###########
###differential RBP expression calculation###
exonsnppair=unique(uniquejoblist[,c("Exon","SNP")])
#filter out SNPs which are not single base 
allele0=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",3)
allele1=sapply(as.character(exonsnppair[,"SNP"]),vectorsplit,"_",4)
allele0length=sapply(allele0,nchar)
allele1length=sapply(allele1,nchar)
singlebase=intersect(which(allele0length==1),which(allele1length==1))
singlebaseexonsnppair=exonsnppair[singlebase,]

sQTLdiffRBP=matrix(NA,dim(singlebaseexonsnppair)[1],dim(edataRBP)[1])
rownames(sQTLdiffRBP)=paste(singlebaseexonsnppair[,1],singlebaseexonsnppair[,2],sep="~")
colnames(sQTLdiffRBP)=rownames(edataRBP)

sQTLdiffRBP_fc=matrix(NA,dim(singlebaseexonsnppair)[1],dim(edataRBP)[1])
rownames(sQTLdiffRBP_fc)=paste(singlebaseexonsnppair[,1],singlebaseexonsnppair[,2],sep="~")
colnames(sQTLdiffRBP_fc)=rownames(edataRBP)

sQTLdiffRBP_dr=matrix(NA,dim(singlebaseexonsnppair)[1],dim(edataRBP)[1])
rownames(sQTLdiffRBP_dr)=paste(singlebaseexonsnppair[,1],singlebaseexonsnppair[,2],sep="~")
colnames(sQTLdiffRBP_dr)=rownames(edataRBP)

#fill in the coefficient and p value matrix
for (i in 1:dim(sQTLdiffRBP)[1]){      
  exonshortID=as.character(singlebaseexonsnppair[i,"Exon"])
  snpID=as.character(singlebaseexonsnppair[i,"SNP"])
  #get the regions that this pair is significant in
  sigregion=as.character(uniquejoblist[rownames(singlebaseexonsnppair)[i],"Brain.region"])
  num.sigregion=length(strsplit(sigregion,split=", ")[[1]])
  
  if (num.sigregion<13){
    for (j in 1:dim(sQTLdiffRBP)[2]){
      rbpID=colnames(sQTLdiffRBP)[j]
      #group the expression TPM by significant/insignificant brain regions
      sigTPM=as.numeric(as.matrix(edataRBP[rbpID,colnames(brainregion)[which(brainregion %in% strsplit(sigregion,split=", ")[[1]])]]))
      insigTPM=as.numeric(as.matrix(edataRBP[rbpID,setdiff(colnames(brainregion),colnames(brainregion)[which(brainregion %in% strsplit(sigregion,split=", ")[[1]])])]))
      difftest=wilcox.test(sigTPM,insigTPM)
      sQTLdiffRBP[i,j]=difftest$p.value
      if (mean(sigTPM,na.rm=T)>mean(insigTPM,na.rm=T)){
        sQTLdiffRBP_fc[i,j]=mean(sigTPM,na.rm=T)/mean(insigTPM,na.rm=T)
        sQTLdiffRBP_dr[i,j]=1
      }else{
        sQTLdiffRBP_fc[i,j]=mean(insigTPM,na.rm=T)/mean(sigTPM,na.rm=T)
        sQTLdiffRBP_dr[i,j]=-1
      }
    }
  }
}

setwd(outputpath)
write.table(sQTLdiffRBP,paste(splicetype,"_",type,"_sQTL_differential_RBP_pvalue.txt",sep=""),sep="\t")
write.table(sQTLdiffRBP_fc,paste(splicetype,"_",type,"_sQTL_differential_RBP_fold_change.txt",sep=""),sep="\t")
write.table(sQTLdiffRBP_dr,paste(splicetype,"_",type,"_sQTL_differential_RBP_direction.txt",sep=""),sep="\t")

sQTLdiffRBP_FDR=matrix(p.adjust(as.numeric(sQTLdiffRBP),"BH"),dim(sQTLdiffRBP)[1],dim(sQTLdiffRBP)[2])
rownames(sQTLdiffRBP_FDR)=rownames(sQTLdiffRBP)
colnames(sQTLdiffRBP_FDR)=colnames(sQTLdiffRBP)
write.table(sQTLdiffRBP_FDR,paste(splicetype,"_",type,"_sQTL_differential_RBP_FDR.txt",sep=""),sep="\t")


###################################################################################################################################################
###########
#section 2#
###########
FDRcutoff=0.05
foldchangecutoff=1

###output significant result###
sigdiffpair=matrix(,nrow=0,ncol=13)
colnames(sigdiffpair)=c("exon.symbol","snp","RBP.symbol","pvalue","FDR","fold.change","direction","sig.region","num.sig.region","exon.fullID","exon.shortID","RBP.fullID","DB.RBPID")

for (i in 1:dim(sQTLdiffRBP)[1]){     
  print(i)
  exonshortID=as.character(singlebaseexonsnppair[i,"Exon"])
  snpID=as.character(singlebaseexonsnppair[i,"SNP"])
  #get the regions that this pair is significant in
  sigregion=as.character(uniquejoblist[rownames(singlebaseexonsnppair)[i],"Brain.region"])
  num.sigregion=length(strsplit(sigregion,split=", ")[[1]])
  
  if (num.sigregion<13){
    for (j in 1:dim(sQTLdiffRBP)[2]){
      rbpID=colnames(sQTLdiffRBP)[j]
      #group the expression TPM by significant/insignificant brain regions
      p.value=sQTLdiffRBP[i,j]
      fold.change=sQTLdiffRBP_fc[i,j]
      direction=sQTLdiffRBP_dr[i,j]
      fdr=sQTLdiffRBP_FDR[i,j]
      dbRBPID=RBPIDconversion[j,"DeepBindID"]

      if (fdr<=FDRcutoff && fold.change>=foldchangecutoff){
        newline=c(as.character(exonIDconversion[exonshortID,"gene.symbol"]),     #exon.symbol
                  snpID,     #snp
                  RBPIDconversion[which(RBPIDconversion[,"edataID"] %in% rbpID),"GeneSymbol"][1],    #RBP.symbol
                  p.value,
                  fdr,
                  fold.change,
                  direction,
                  sigregion,
                  num.sigregion,
                  as.character(exonIDconversion[exonshortID,"fullID"]),   #exon.fullID
                  as.character(exonIDconversion[exonshortID,"shortID"]),  #exon.shortID
                  rbpID,      #RBP full ID
                  dbRBPID)
        sigdiffpair=rbind(sigdiffpair,newline)
      }
    }
  }
}
setwd(outputpath)
write.table(sigdiffpair,paste(splicetype,"_",type,"_significant_sQTL_RBP_pairs.txt",sep=""),sep="\t",row.names=F)


