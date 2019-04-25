expressiontype="allgene"

print(expressiontype)

library(boot)
#########################
#read in gene expression#
#########################
inputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_expression_TPM"
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
IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

print(dim(expr))

####################
#read in annotation#
####################
phenotypepath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"
setwd(phenotypepath)
oriage=read.table("agetable_brain.txt",sep="\t",header=T)           
origender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
tempgender=origender
origender[1,which(as.numeric(as.matrix(tempgender)[1,])==1)]="male"
origender[1,which(as.numeric(as.matrix(tempgender)[1,])==2)]="female"
oribrainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)

#get the phenotype information of samples in the expr matrix (not all 1409 sampes are in the expr data)
age=oriage[,as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"])]
gender=origender[,as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"])]
brainregion=oribrainregion[,as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"])]

print(dim(age))

#####
#SVA#
#####
pheno=cbind(t(as.matrix(gender)),t(as.matrix(brainregion)),t(as.matrix(age)))   #phenotypes on the column and samples on the row
pheno=as.data.frame(pheno)
edata = as.matrix(expr[,-c(1,2)])
edata=edata[which(apply(t(edata), 2, var)!=0),]   #remove genes with 0 variance
#log transformation of TPM
edata=log2(edata+10^-10)

outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/8.2_gene_expression_SVA/result/",expressiontype,sep="")
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)

library(sva)

###brain region + gender + age###
mod = model.matrix(~as.factor(brain_region) + as.factor(age) + as.factor(gender), data=pheno)
mod0 = model.matrix(~1,data=pheno)

n.sv = num.sv(edata,mod,method="leek")
n.sv

svobj = sva(edata,mod,mod0,n.sv=n.sv)

if(n.sv>0){
  svobj$sv=as.matrix(svobj$sv)
  rownames(svobj$sv)=colnames(age)
  write.table(svobj$sv,paste(expressiontype,"_SVA_brain_gender_age_surrogate_variable.txt",sep=""),sep="\t")
  write.table(svobj$pprob.gam,paste(expressiontype,"_SVA_brain_gender_age_pprob.gam.txt",sep=""),sep="\t")
  write.table(svobj$pprob.b,paste(expressiontype,"_SVA_brain_gender_age_pprob.b.txt",sep=""),sep="\t")
}





