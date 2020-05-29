splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

library(boot)

inputpath="./02_PSI_value_quantification/02_Imputation/example_output"
setwd(inputpath)
impute_PSI_filter=read.table(paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t",header=T)
#transform PSI to logit(PSI)
transformedpsi=as.numeric(as.matrix(impute_PSI_filter))
transformedpsi[which(transformedpsi<0.01)]=0.01
transformedpsi[which(transformedpsi>0.99)]=0.99
transformedpsi=logit(transformedpsi)
logitPSI=matrix(transformedpsi,dim(impute_PSI_filter)[1],dim(impute_PSI_filter)[2])
rownames(logitPSI)=rownames(impute_PSI_filter)
colnames(logitPSI)=colnames(impute_PSI_filter)
logitPSI=as.data.frame(logitPSI)
impute_PSI_filter=logitPSI

#read in phenotype information
phenotypepath="./03_Get_sample_annotation/example_output"
setwd(phenotypepath)
age=read.table("agetable_brain.txt",sep="\t",header=T)           
gender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
tempgender=gender
gender[1,which(as.numeric(as.matrix(tempgender)[1,])==1)]="male"
gender[1,which(as.numeric(as.matrix(tempgender)[1,])==2)]="female"
brainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)

pheno=cbind(t(as.matrix(gender)),t(as.matrix(brainregion)))   #phenotypes on the column and samples on the row
pheno=as.data.frame(pheno)
edata = as.matrix(impute_PSI_filter)

outputpath="/output/path"
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)



#####
#SVA#
#####
library(sva)

###only brain region###
mod = model.matrix(~as.factor(brain_region), data=pheno)
mod0 = model.matrix(~1,data=pheno)

#method 1 (this is based on the SVA document)
n.sv = num.sv(edata,mod,method="leek")
n.sv

svobj = sva(edata,mod,mod0,n.sv=n.sv)

if(n.sv>0){
  svobj$sv=as.matrix(svobj$sv)
  rownames(svobj$sv)=colnames(age)
  write.table(svobj$sv,paste(splicetype,"_",counttype,"_SVA_brain_surrogate_variable.txt",sep=""),sep="\t")
  write.table(svobj$pprob.gam,paste(splicetype,"_",counttype,"_SVA_brain_pprob.gam.txt",sep=""),sep="\t")
  write.table(svobj$pprob.b,paste(splicetype,"_",counttype,"_SVA_brain_pprob.b.txt",sep=""),sep="\t")
}



