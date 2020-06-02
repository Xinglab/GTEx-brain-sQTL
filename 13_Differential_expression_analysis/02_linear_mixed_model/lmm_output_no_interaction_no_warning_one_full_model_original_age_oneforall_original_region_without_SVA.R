expressiontype="allgene"

library(boot)
library(lme4)  
library(ResourceSelection)
#########################
#read in gene expression#
#########################
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
edata = as.matrix(expr[,-c(1,2)])

IDconversion=read.table("SRRID_2_GTEXID.txt",sep="\t",header=T)
rownames(IDconversion)=IDconversion[,"sampleID"]

####################
#read in annotation#
####################
phenotypepath="./03_Get_sample_annotation/example_output"
setwd(phenotypepath)
orisampleID_SRRID=read.table("sampleID_SRRID_brain.txt",sep="\t",header=T)
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
sampleID_SRRID=orisampleID_SRRID[as.character(IDconversion[colnames(expr)[-c(1,2)],"SRRlist"]),]

####################
#read in SVA result#
####################
numSV=0             #we don't correct for SVs as explained in the manuscript

###########
#functions#
###########
lmm.LRT <- function(data) {
  options(warn=2) #disable the print of error and warning information. 
  #If we use this, models with warnings will be considered as not converged.
  #If we don't use this, as long as there is no error, the model will be considered as converged, even though sometimes the warning says the model failed to converge. 
  
  age.pval <- NA     #p value=NA if model doesn't converge
  gender.pval <- NA  #p value=NA if model doesn't converge
  br.pval <- NA  #p value=NA if model doesn't converge
  betas <-  rep(0,(19+data$nsv))      #coefficients of fixed effect if model doesn't converge
  
  response <- data$x
  char=as.character(sampleID_SRRID[,"sampleID"])   #SRRID and sampleID is in the same order with all the other input variables (age, gender, etc)
  num=seq(1:length(unique(char)))    #each unique sample ID will have a number
  temp=cbind(unique(char),num)       #we assign each number to a donor
  rownames(temp)=unique(char)
  colnames(temp)=c("char","num")
  obs <- as.numeric(temp[char,"num"])    #obs can be either numeric or factor   
  
  if (data$nsv==0){
    testglm = try(suppressMessages(lmer(response ~ data$AGE + data$GENDER + data$BR + (1|obs) , REML=FALSE)),silent=TRUE)  #invidual-level random effect for overdispersion. 
    if (!(inherits(testglm,"try-error"))) { # if the full model converges. if converge, inherits(XXX,"try-error")=FALSE, otherwise TRUE. 
      testglmAGE = try(suppressMessages(lmer(response ~ data$GENDER + data$BR + (1|obs), REML=FALSE)),silent=TRUE)
      testglmGENDER= try(suppressMessages(lmer(response ~ data$AGE + data$BR + (1|obs), REML=FALSE)),silent=TRUE)
      testglmBR=try(suppressMessages(lmer(response ~ data$AGE + data$GENDER + (1|obs), REML=FALSE)),silent=TRUE)
      
      if (length(fixef(testglm))==length(betas)){ 
         betas <- fixef(testglm)
      }
      
      if (!(inherits(testglmAGE,"try-error"))){   #if the age model also converges
        options(error = NULL, warn = 0)       #there may be some warnings here and we don't treat them as error
        age.pval <-  anova(testglm, testglmAGE)$"Pr(>Chisq)"[2]
        options(warn=2)
      }
      if (!(inherits(testglmGENDER,"try-error"))){   #if the gender model also converges
        options(error = NULL, warn = 0)
        gender.pval <-  anova(testglm, testglmGENDER)$"Pr(>Chisq)"[2]
        options(warn=2)
      }
      if (!(inherits(testglmBR,"try-error"))){   #if the brain region model also converges
        options(error = NULL, warn = 0)
        br.pval <-  anova(testglm, testglmBR)$"Pr(>Chisq)"[2]
        options(warn=2)
      }
    }
  }

  options(error = NULL, warn = 0)
  return(list(betas=betas, pval=c(age.pval,gender.pval,br.pval)))  
}


betamatrix=matrix(NA,dim(edata)[1],(19+numSV))      
rownames(betamatrix)=rownames(edata)

pvmatrix=matrix(NA,dim(edata)[1],3)
rownames(pvmatrix)=rownames(edata)
colnames(pvmatrix)=c("pAGE","pGENDER","pBR")

for (gene in 1:dim(edata)[1]){
  #log transformation of TPM
  EXPR=log2(edata[gene,]+10^-10)
  
  if (numSV==0){
    lmminput=list(x=EXPR,AGE=as.character(as.matrix(age[1,])),GENDER=as.character(gender[1,]),BR=as.character(as.matrix(brainregion[1,])),nsv=numSV)
  }

  test=lmm.LRT(lmminput)
  colnames(betamatrix)=names(test$betas)
  betamatrix[gene,]=test$betas
  pvmatrix[gene,]=test$pval
}


########
#output#
########
outputpath="/output/path"
label=strsplit(outputpath,split="/")[[1]][12]
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)
write.table(betamatrix,paste("betamatrix_",label,".txt",sep=""),sep="\t",quote=F)
write.table(pvmatrix,paste("pvmatrix_",label,".txt",sep=""),sep="\t",quote=F)
