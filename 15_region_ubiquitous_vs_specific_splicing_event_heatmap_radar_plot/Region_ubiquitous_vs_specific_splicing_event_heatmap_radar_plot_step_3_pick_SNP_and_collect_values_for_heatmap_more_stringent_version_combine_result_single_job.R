rangelist=c("200kb")
filelist=c("pvalue","beta")
outputprefixlist=c("pvmatrix","betamatrix")
rootinput="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_3_pick_SNP_and_collect_values_for_heatmap_more_stringent_version/single_job_run/result/200kb"
sqtlrootinput=paste("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/sQTL_run/logit/JC",
                    "SE",sep="/")

###################################
#summarize exon result into matrix#
###################################
for (windowsize in rangelist){
  for (f in 1:length(filelist)){
    file=filelist[f]
    
    inputpath=paste(rootinput,"/",file,sep="")
    
    #merge the file for each exon into a big table
    setwd(inputpath)
    command=paste("cat *.txt > ",rootinput,"/",paste(outputprefixlist[f],"_",windowsize,".txt",sep=""),sep="")
    system(command)
  }
}

############################
#add column names to matrix#
############################
brainregionlist=c("Brain-Amygdala",
                  "Brain-AnteriorcingulatecortexBA24",
                  "Brain-Caudatebasalganglia",
                  "Brain-CerebellarHemisphere",
                  "Brain-Cerebellum",
                  "Brain-Cortex",
                  "Brain-FrontalCortexBA9",
                  "Brain-Hippocampus",
                  "Brain-Hypothalamus",
                  "Brain-Nucleusaccumbensbasalganglia",
                  "Brain-Putamenbasalganglia",
                  "Brain-Spinalcordcervicalc-1",
                  "Brain-Substantianigra")


#200kb
setwd(rootinput)
betamatrix200k=read.table("betamatrix_200kb.txt",sep="\t",row.names = 1)
pvmatrix200k=read.table("pvmatrix_200kb.txt",sep="\t",row.names = 1)

colnames(betamatrix200k)=colnames(pvmatrix200k)=brainregionlist

write.table(betamatrix200k,"betamatrix_200kb.txt",sep="\t")
write.table(pvmatrix200k,"pvmatrix_200kb.txt",sep="\t")


############################################################################################
#calculate number of significant regions & the location of each SNP & the rs ID of each SNP#
############################################################################################
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/200kb_highp")
lcmatrix200k=read.table("lcmatrix200k_highp.txt",sep="\t",header=T,check.names=F)
topsnpmatrix200k=read.table("topsnpmatrix200k_highp.txt",sep="\t",header=T,check.names=F)
rsmatrix200k=read.table("rsmatrix200k_highp.txt",sep="\t",header=T,check.names=F)

cutoff.p=1e-5

#get the permutation p value cutoff for each brain region
cutoff.fdr=matrix(NA,length(brainregionlist),1)
rownames(cutoff.fdr)=brainregionlist
for (br in 1:length(brainregionlist)){
  currentbr=brainregionlist[br]
  setwd(paste(sqtlrootinput,"/",currentbr,sep=""))
  cutoff.fdr[br,1]=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
}

numsigregion.p=numsigregion.fdr=lclist=rslist=rep(NA,dim(pvmatrix200k)[1])
for (i in 1:dim(pvmatrix200k)[1]){
  numsigregion.p[i]=sum(pvmatrix200k[i,]<cutoff.p)
  numsigregion.fdr[i]=sum(pvmatrix200k[i,]<cutoff.fdr)
  temp=rownames(pvmatrix200k)[i]
  currentexon=strsplit(temp,split="~")[[1]][1]
  currentsnp=strsplit(temp,split="~")[[1]][2]
  lclist[i]=as.character(as.matrix(lcmatrix200k[currentexon,which(topsnpmatrix200k[currentexon,]==currentsnp)[1]]))
  rslist[i]=as.character(as.matrix(rsmatrix200k[currentexon,which(topsnpmatrix200k[currentexon,]==currentsnp)[1]]))
}

setwd(rootinput)
output=cbind(numsigregion.p,numsigregion.fdr,lclist,rslist)
rownames(output)=rownames(pvmatrix200k)
write.table(output,"numsigregion_lc_rs_list200k.txt",sep="\t")









