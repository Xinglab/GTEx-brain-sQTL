#The purpose of this code is to get the list of age, gender, and brain region dependent genes based on FDR<5%

expressiontype="allgene"

label="lmm_output_no_interaction_no_warning_one_full_model_original_age_oneforall_original_region_without_SVA"
rootinput="/output/path/of/lmm/result"
fdrcutoff=0.05

###################
#read in phenotype#
###################
inputpath="/input/to/GTEx/processed/gene/expression/data"
setwd(inputpath)
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


###########################
#get the significant genes#
###########################
setwd(rootinput)
pvfilename=paste("pvmatrix_",label,".txt",sep="")
pvmatrix=read.table(pvfilename,sep="\t",header=T)
rownames(pvmatrix)=rownames(edata)

#get the events information from the row name
eventsinfo=matrix(NA,dim(pvmatrix)[1],2)
for (i in 1:dim(pvmatrix)[1]){
  temp=strsplit(rownames(pvmatrix)[i],split="_")[[1]]
  eventsinfo[i,]=temp[1:2]
}

#calculate FDR
temp=as.numeric(as.matrix(pvmatrix))
temp[which(is.na(temp))]=1           #assign 1 to missing p values (not converged)
fdrmatrix=matrix(temp,dim(pvmatrix)[1],dim(pvmatrix)[2])
rownames(fdrmatrix)=rownames(pvmatrix)
colnames(fdrmatrix)=colnames(pvmatrix)
for (i in 1:dim(fdrmatrix)[2]){
  fdrmatrix[,i]=p.adjust(fdrmatrix[,i],"BH")
}

#get events to plot
fdrDS.BR=which(fdrmatrix[,"pBR"]<fdrcutoff)    #events pass FDR cutoff
fdrDS.AGE=which(fdrmatrix[,"pAGE"]<fdrcutoff) 
fdrDS.GENDER=which(fdrmatrix[,"pGENDER"]<fdrcutoff) 

setwd("/output/path")
if (length(fdrDS.BR)>0){
  outputtable=cbind(eventsinfo[fdrDS.BR,],pvmatrix[fdrDS.BR,"pBR"],fdrmatrix[fdrDS.BR,"pBR"])
  colnames(outputtable)=c("Gene ID","Gene Symbol",
                          "P-value","FDR")
  write.table(outputtable,paste("BR_dependent_gene_summary.txt",sep=""),sep="\t")
}
if (length(fdrDS.AGE)>0){
  outputtable=cbind(eventsinfo[fdrDS.AGE,],pvmatrix[fdrDS.AGE,"pAGE"],fdrmatrix[fdrDS.AGE,"pAGE"])
  colnames(outputtable)=c("Gene ID","Gene Symbol",
                          "P-value","FDR")
  write.table(outputtable,paste("AGE_dependent_gene_summary.txt",sep=""),sep="\t")
}
if (length(fdrDS.GENDER)>0){
  outputtable=cbind(eventsinfo[fdrDS.GENDER,],pvmatrix[fdrDS.GENDER,"pGENDER"],fdrmatrix[fdrDS.GENDER,"pGENDER"])
  colnames(outputtable)=c("Gene ID","Gene Symbol",
                          "P-value","FDR")
  write.table(outputtable,paste("GENDER_dependent_gene_summary.txt",sep=""),sep="\t")
}



