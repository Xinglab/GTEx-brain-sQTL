#The purpose of this code is just to plot the correlation between Expression ratio (T-STAR:Sam68) and PSI (NRXN2)
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
#edata=edata[which(apply(t(edata), 2, var)!=0),]   #remove genes with 0 variance
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


##############################
#get the RBP expression ratio#
##############################
geneexprratio=edata["ENSG00000131773.9_KHDRBS3",]/edata["ENSG00000121774.13_KHDRBS1",]   
#ENSG00000121774.13_KHDRBS1      #Sam68
#ENSG00000131773.9_KHDRBS3       #T-STAR
 
###################
#get the PSI value#
###################
exonpsi=PSI["30935|ENSG00000110076.18_2|NRXN2|chr11|-|64393934|64394024|64390224|64390550|64397873|64398045",]


##################
#correlation plot#
##################
library(LSD)
setwd(outputpath)
pdf("correlation between ratio and NRXN2 psi.pdf",height=6,width=6)
heatscatter(as.numeric(geneexprratio), as.numeric(exonpsi), xlab='Expression ratio (T-STAR:Sam68)', ylab='PSI (NRXN2)', add.contour=F)
dev.off()


