splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)

inputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/4.2_data_inspection_imputation/result/",counttype,sep="")
setwd(inputpath)

impute_PSI_filter=read.table(paste(splicetype,"_",counttype,"_impute_PSI_filter_KNN.txt",sep=""),sep="\t",header=T)

###logit transformation###
library(boot)
transformedpsi=as.numeric(as.matrix(impute_PSI_filter))
transformedpsi[which(transformedpsi<0.01)]=0.01
transformedpsi[which(transformedpsi>0.99)]=0.99
transformedpsi=logit(transformedpsi)
logitPSI=matrix(transformedpsi,dim(impute_PSI_filter)[1],dim(impute_PSI_filter)[2])
rownames(logitPSI)=rownames(impute_PSI_filter)
colnames(logitPSI)=colnames(impute_PSI_filter)
logitPSI=as.data.frame(logitPSI)
impute_PSI_filter=logitPSI
###logit transformation###

outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/",counttype,sep="")
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)


library(Rtsne)

###t-sne on all exons###
#clustering of samples (the clustering of t-SNE is on the rows of a given matrix and the row of the transposed matrix is sample)
tsne_out <- Rtsne(t(as.matrix(impute_PSI_filter)))
temp=tsne_out$Y
rownames(temp)=rownames(t(impute_PSI_filter))

write.table(temp,paste(splicetype,"_",counttype,"_tsne_sample_logit_PSI.txt",sep=""),sep="\t")      ###logit PSI###





