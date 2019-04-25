splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

inputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/1_get_matrix_from_rMATS/result/",counttype,sep="")
setwd(inputpath)

PSI_filter=read.table(paste(splicetype,"_",counttype,"_PSI_filter.txt",sep=""),sep="\t",header=T)
#impute missing value (usually we don't like impute data but the imputed data is only used in PCA. We will use the original data in our GLMM model)
library(impute)
imputation=impute.knn(as.matrix(PSI_filter), k=30, maxp=dim(PSI_filter)[1])
impute_PSI_filter <- imputation$data

impute_PSI_filter[which(impute_PSI_filter>1)]=1
impute_PSI_filter[which(impute_PSI_filter<0)]=0

outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.2_data_inspection_imputation/result/",counttype,sep="")
command=paste("mkdir -p ",outputpath,sep="")
system(command)

setwd(outputpath)
write.table(impute_PSI_filter,paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t")

