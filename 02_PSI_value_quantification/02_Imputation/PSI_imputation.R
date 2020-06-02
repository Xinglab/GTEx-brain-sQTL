splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

inputpath="./02_PSI_value_quantification/01_Get_PSI_from_rMATS_output/example_output"
setwd(inputpath)

PSI_filter=read.table(paste(splicetype,"_",counttype,"_PSI_filter.txt",sep=""),sep="\t",header=T)

library(impute)
imputation=impute.knn(as.matrix(PSI_filter), k=30, maxp=dim(PSI_filter)[1])
impute_PSI_filter <- imputation$data

impute_PSI_filter[which(impute_PSI_filter>1)]=1
impute_PSI_filter[which(impute_PSI_filter<0)]=0

outputpath="./02_PSI_value_quantification/02_Imputation/example_output"
command=paste("mkdir -p ",outputpath,sep="")
system(command)

setwd(outputpath)
write.table(impute_PSI_filter,paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t")

