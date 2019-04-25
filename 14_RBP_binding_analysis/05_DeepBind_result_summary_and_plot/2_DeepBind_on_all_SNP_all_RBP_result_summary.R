splicetypelist=c("SE","A3SS","A5SS")
typelist=c("pvalue","permutation")

for (splicetype in splicetypelist){
  for (type in typelist){
    #generate output folder
    inputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_plot/",splicetype,"/",type,sep="")
    outputpath=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_summary/",splicetype,"/",type,sep="")
    command=paste("mkdir -p",outputpath)
    system(command)
    
    #get all the folders under the current folder
    setwd(inputpath)
    folderlist=system("ls",intern=T)
    for (f in 1:length(folderlist)){
      #in each folder, we have the result for all RBPs, we merge them into one file first and output it to the output folder
      command=paste("cat ./",folderlist[f],"/*/summarystats_*.txt > ",outputpath,"/summary_",f,".txt",sep="")
      system(command)
    }
    
    #then we merge all the sub files together
    setwd(outputpath)
    command=paste("cat summary_*.txt > ",outputpath,"/all_summary.txt",sep="")
    system(command)
    
    
    #looking for missing files
    #read in all the jobs in the summary file
    setwd(outputpath)
    result=read.table("all_summary.txt",sep="\t")
    rownames(result)=result[,1]
    result=result[,-1]
    colnames(result)=c("single.SNP.value","average.SNP.value","single.SNP.rank","single.base.height","average.base.height","single.base.rank","outputpath")
    
    #read in the original job list
    rootoutput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
    setwd(rootoutput)
    uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")
    
    missing=uniquejoblist[setdiff(rownames(uniquejoblist),rownames(result)),]
    setwd(outputpath)
    write.table(missing,"missing_summary.txt",sep="\t")
  }
}


