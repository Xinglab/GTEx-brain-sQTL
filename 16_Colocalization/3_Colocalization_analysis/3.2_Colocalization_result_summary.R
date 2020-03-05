#the purpose of this code is to summarize all the colocalization result

inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis/result"
outputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis"

setwd(inputpath)
command=paste("cat *.txt > ",outputpath,"/colocalization_result_summary.txt",sep="")
system(command)

setwd(outputpath)
result=read.table("colocalization_result_summary.txt",sep="\t")

rownames(result)=result[,1]
result=result[,-1]

colnames(result)=c("sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF","sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF","sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF",
                   "sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF","sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF","sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF")

setwd(outputpath)
write.table(result,"colocalization_result_summary.txt",sep="\t")
  
  
#####################################################################################################################
inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis"
setwd(inputpath)
result=read.table("colocalization_result_summary.txt",sep="\t",header=T)

for (i in 1:6){
  for (j in c(0.75,0.8)){
    print(paste(sum(result[,i]>=j,na.rm=T),"events with PP4 greater or equal to",j,"based on",colnames(result)[i]))
  }
}

##############################################
#requiring_sample_size_in_summary_statistics:#
##############################################
#[1] "14 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"
#[1] "11 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"

#[1] "14 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"
#[1] "11 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"

#[1] "18 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"
#[1] "15 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"

##################################################
#not requiring_sample_size_in_summary_statistics:#
##################################################
#[1] "64 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"
#[1] "53 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"

#[1] "64 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"
#[1] "53 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"

#[1] "87 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"
#[1] "72 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"

#[1] "22 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF"
#[1] "19 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF"

#[1] "22 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF"
#[1] "19 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF"

#[1] "33 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF"
#[1] "29 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF"



