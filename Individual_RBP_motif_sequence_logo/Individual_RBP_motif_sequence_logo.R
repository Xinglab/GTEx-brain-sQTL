library(ggseqlogo)
library(ggplot2)
library(reticulate)
functionfile="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/motif_scan/rbp_map_functions.py"
source_python(functionfile)

outputpath="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/10.6_individual_RBP_motif_sequence_logo/result"

#read in the motif information
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/data/motif_information")
motiflist=read.table("RBP_Binding_Sites.txt",sep="\t",header=T,skipNul=TRUE)

RBPlist=c("RBFOX1","NOVA1","HNRNPK")

for (i in 1:length(RBPlist)){
  RBP=RBPlist[i]
  RBPmotif=as.character(as.matrix(motiflist)[which(motiflist[,1] %in% RBP),"Consensus.binding.site"])
  mlist=return_seq_log_list(RBPmotif)[[1]]
  
  #generate the plot
  seq2plot=mlist
  setwd(outputpath)
  outfile=paste("Seqlogo_",RBP,".pdf",sep="")
  pdf(outfile,width=4,height=4)
  p=ggseqlogo(seq2plot,method = 'prob')
  print(p)
  dev.off()
}


