expressiontype="allgene"  

#########################
#read in gene expression#
#########################
inputpath="/input/to/GTEx/processed/gene/expression/data"
setwd(inputpath)
if (expressiontype=="allgene"){
  expr=read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_brain_normalized.txt",header=T,sep="\t",check.names=F)     #normalized TPM
}
if (expressiontype=="allRBP"){
  expr=read.table("all_RBP_normalized_expression.txt",header=T,sep="\t",check.names=F)     #normalized TPM
}
if (expressiontype=="allSF"){
  expr=read.table("SF_normalized_expression.txt",header=T,sep="\t",check.names=F)     #normalized TPM
}

print(dim(expr))

edata = as.matrix(expr[,-c(1,2)])
edata=edata[which(apply(t(edata), 2, var)!=0),]   #remove genes with 0 variance
#log transformation of TPM
edata=log2(edata+10^-10)        



outputpath="/output/path"
command=paste("mkdir -p ",outputpath,sep="")
system(command)
setwd(outputpath)


library(Rtsne)

tsne_out <- Rtsne(t(as.matrix(edata)))
temp=tsne_out$Y
rownames(temp)=rownames(t(edata))

write.table(temp,paste(expressiontype,"_tsne_sample_log_expr.txt",sep=""),sep="\t")   





