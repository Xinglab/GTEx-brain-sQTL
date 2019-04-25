code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/7.3_DS_analysis_oneforalll_linear_mixed_model/lmm_output_no_interaction_no_warning_one_full_model_original_age_imputed_psi_no_weight_JC_oneforall_original_region_with_SVA"
splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

#######
#input#
#######
library(boot)
library(lme4)  
library(ResourceSelection)

inputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.2_data_inspection_imputation/result/",counttype,sep="")
setwd(inputpath)
impute_PSI_filter=read.table(paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t",header=T)
psi=impute_PSI_filter

#read in annotation information
phenotypepath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/2_get_phenotype_information/result"
setwd(phenotypepath)
age=read.table("agetable_brain.txt",sep="\t",header=T)           
gender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
brainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)
sampleID_SRRID=read.table("sampleID_SRRID_brain.txt",sep="\t",header=T)

###read in SVA result###
inputSVApath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/7.1_DS_analysis_SVA/result/",counttype,sep="")
inputSVAname=paste(splicetype,"_",counttype,"_SVA_brain_gender_age_surrogate_variable.txt",sep="")
setwd(inputSVApath)
sva=try(suppressMessages(read.table(inputSVAname,sep="\t")),silent=TRUE)      #if there is no SVA result, it will generate a warning. Don't worry about that. It won't impact the result. 
if (inherits(sva,"try-error")){       
  numSV=0    #if we don't have this file, it means that there is no surrogate variables
}else{
  numSV=dim(sva)[2]
}

print(numSV)

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
  #obs <- seq(1,length(data$y))    #if each sample is from a differnet persion, then we can use this 
  #if there are samples from the same person, we need to do it as following:
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

  if (data$nsv==2){
    testglm = try(suppressMessages(lmer(response ~ data$AGE + data$GENDER + data$BR + data$SV1 + data$SV2 + (1|obs) , REML=FALSE)),silent=TRUE)  #invidual-level random effect for overdispersion. 
    if (!(inherits(testglm,"try-error"))) { # if the full model converges. if converge, inherits(XXX,"try-error")=FALSE, otherwise TRUE. 
      testglmAGE = try(suppressMessages(lmer(response ~ data$GENDER + data$BR + data$SV1 + data$SV2 + (1|obs), REML=FALSE)),silent=TRUE)
      testglmGENDER= try(suppressMessages(lmer(response ~ data$AGE + data$BR + data$SV1 + data$SV2 + (1|obs), REML=FALSE)),silent=TRUE)
      testglmBR=try(suppressMessages(lmer(response ~ data$AGE + data$GENDER + data$SV1 + data$SV2 + (1|obs), REML=FALSE)),silent=TRUE)
      
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



betamatrix=matrix(NA,dim(psi)[1],(19+numSV))      
rownames(betamatrix)=rownames(psi)

pvmatrix=matrix(NA,dim(psi)[1],3)
rownames(pvmatrix)=rownames(psi)
colnames(pvmatrix)=c("pAGE","pGENDER","pBR")

for (exon in 1:dim(psi)[1]){
  #because we are doing linear mixed model, we need to use logit(psi). logit(0) and logit(1) will generate Inf/-Inf so we truncate the original psi
  transformedpsi=as.numeric(psi[exon,])
  transformedpsi[which(transformedpsi<0.01)]=0.01     #there may be missing values but it will not influence the result
  transformedpsi[which(transformedpsi>0.99)]=0.99
  transformedpsi=logit(transformedpsi)

  if (numSV==0){
    lmminput=list(x=transformedpsi,AGE=as.character(as.matrix(age[1,])),GENDER=as.character(gender[1,]),BR=as.character(as.matrix(brainregion[1,])),nsv=numSV)
  }
  if (numSV==2){
    lmminput=list(x=transformedpsi,AGE=as.character(as.matrix(age[1,])),GENDER=as.character(gender[1,]),BR=as.character(as.matrix(brainregion[1,])),SV1=sva[,1],SV2=sva[,2],nsv=numSV)
  }
  
  test=lmm.LRT(lmminput)
  colnames(betamatrix)=names(test$betas)
  betamatrix[exon,]=test$betas
  pvmatrix[exon,]=test$pval
}


########
#output#
########
outputpath=paste(code_folder,"result",splicetype,counttype,sep="/")
label=strsplit(outputpath,split="/")[[1]][12]
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)
write.table(betamatrix,paste("betamatrix_",label,".txt",sep=""),sep="\t",quote=F)
write.table(pvmatrix,paste("pvmatrix_",label,".txt",sep=""),sep="\t",quote=F)
