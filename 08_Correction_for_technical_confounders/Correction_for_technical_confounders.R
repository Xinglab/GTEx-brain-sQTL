splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

###read in original PSI value and transform that into logit PSI###
library(boot)
inputpath="./02_PSI_value_quantification/02_Imputation/example_output"
setwd(inputpath)
impute_PSI_filter=read.table(paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t",header=T)   #we use imputed PSI here since the SVs are based on SVA analysis and the SVA analysis is done on imputed PSI value (SVA cannot be run on data with missing values)
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

outputpath="/output/path"
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)


###read in genotype PCA result###
setwd("/path/for/genotype/PCA/result")
inputvec="PCA_result_prefix.eigenvec"
pca.x=read.table(inputvec,sep=" ")
rownames(pca.x)=pca.x[,1]
pca.x=pca.x[,-c(1,2)]
colnames(pca.x)=paste("PC",seq(1,dim(pca.x)[2]),sep="")

#read in ID conversion information (PCA is done in donors but all the other results are based on samples)
sample.pca.x=matrix(NA,dim(impute_PSI_filter)[2],dim(pca.x)[2])
rownames(sample.pca.x)=colnames(impute_PSI_filter)
colnames(sample.pca.x)=colnames(pca.x)

IDconversionpath="/input/path/to/GTEx/brain/tissue/sample/annotation"
IDconversionname="gtex_v7_brain.csv" 
setwd(IDconversionpath)
IDconversion=read.csv(IDconversionname,header=T)

for (i in 1:dim(sample.pca.x)[1]){
  SRRID=rownames(sample.pca.x)[i]
  sampleID=as.character(IDconversion[which(IDconversion[,"Run_s"]==SRRID),"submitted_subject_id_s"])[1]   #there are duplications in the annotation file
  sample.pca.x[i,]=as.matrix(pca.x[sampleID,])   #some samples come from donors that are not in the genotype data because originally the genotype had 652 individuals but some of them (17 individuals) are filtered out at QC step. 
}


###read in SVA result###
inputSVApath="/path/for/SVA/result/on/PSI"
inputSVAname=paste(splicetype,"_",counttype,"_SVA_brain_surrogate_variable.txt",sep="")
setwd(inputSVApath)
SVA=try(suppressMessages(read.table(inputSVAname,sep="\t")),silent=TRUE)
if (inherits(SVA,"try-error")){       
  numSV=0    #if we don't have this file, it means that there is no surrogate variables
}else{
  numSV=dim(SVA)[2]
}


###get new logit PSI value after regressing out confounders###
new_logit_PSI=matrix(NA,dim(impute_PSI_filter)[1],dim(impute_PSI_filter)[2])
rownames(new_logit_PSI)=rownames(impute_PSI_filter)
colnames(new_logit_PSI)=colnames(impute_PSI_filter)

cutoff=3   #we choose the top 3 PC

for (i in 1:dim(new_logit_PSI)[1]){    #for each exon
  y=as.numeric(as.matrix(impute_PSI_filter[i,]))
  if (numSV==0){
    model=lm(y~sample.pca.x[,"PC1"]+sample.pca.x[,"PC2"]+sample.pca.x[,"PC3"])     #we choose the top 3 PC
  }
  if (numSV==1){
    model=lm(y~sample.pca.x[,"PC1"]+sample.pca.x[,"PC2"]+sample.pca.x[,"PC3"]+SVA[,1])
  }
  if (numSV==2){
    model=lm(y~sample.pca.x[,"PC1"]+sample.pca.x[,"PC2"]+sample.pca.x[,"PC3"]+SVA[,1]+SVA[,2])
  }
  if (numSV==3){
    model=lm(y~sample.pca.x[,"PC1"]+sample.pca.x[,"PC2"]+sample.pca.x[,"PC3"]+SVA[,1]+SVA[,2]+SVA[,3])
  }
  new_logit_PSI[i,which(!is.na(sample.pca.x[,"PC1"]))]=resid(model)      #there are missing values in PCA of genotype 
  
  #for nunSV>0, we can simply use model=lm(y~as.matrix(sample.pca.x[,1:3])+as.matrix(SVA)) to represent all cases
  if (numSV>0){
    test=lm(y~as.matrix(sample.pca.x[,1:cutoff])+as.matrix(SVA))
    if(!identical(resid(model),resid(test))){
      print("residual not identical")
    }
  }
}

setwd(outputpath)
write.table(new_logit_PSI,paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep=""),sep="\t")

