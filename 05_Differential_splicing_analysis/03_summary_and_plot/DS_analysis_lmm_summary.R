#The purpose of this code is to get the list of age, gender, and brain region dependent splicing events for SE, A3SS and A5SS based on FDR<5%

splicetypelist=c("SE","A3SS","A5SS")
counttype="JC"
rootinput="./05_Differential_splicing_analysis/02_linear_mixed_model"
fdrcutoff=0.05
outputpath="./05_Differential_splicing_analysis/03_summary_and_plot/example_output"

#read in phenotype information
phenotypepath="./03_Get_sample_annotation/example_output"
setwd(phenotypepath)
oriage=read.table("agetable_brain.txt",sep="\t",header=T)           
gender=read.table("gendertable_brain.txt",sep="\t",header=T)          #current gender is numeric data, not categorical data. It needs to be changed to categorical data for further model fitting
brainregion=read.table("brain_region_table_brain.txt",sep="\t",header=T)
sampleID_SRRID=read.table("sampleID_SRRID_brain.txt",sep="\t",header=T)

agelabel=c("","merged")  #whether we want to plot using original age or merged age (the order must be like this)

#change brain region labels
brainregion[which(brainregion=="Brain - Cerebellar Hemisphere")]="Cerebellar Hemisphere"
brainregion[which(brainregion=="Brain - Cerebellum")]="Cerebellum"
brainregion[which(brainregion=="Brain - Cortex")]="Cortex"
brainregion[which(brainregion=="Brain - Frontal Cortex (BA9)")]="Frontal Cortex"
brainregion[which(brainregion=="Brain - Anterior cingulate cortex (BA24)")]="Anterior cingulate cortex"
brainregion[which(brainregion=="Brain - Hypothalamus")]="Hypothalamus"
brainregion[which(brainregion=="Brain - Caudate (basal ganglia)")]="Caudate"
brainregion[which(brainregion=="Brain - Amygdala")]="Amygdala"
brainregion[which(brainregion=="Brain - Spinal cord (cervical c-1)")]="Spinal cord"
brainregion[which(brainregion=="Brain - Nucleus accumbens (basal ganglia)")]="Nucleus accumbens"
brainregion[which(brainregion=="Brain - Substantia nigra")]="Substantia nigra"
brainregion[which(brainregion=="Brain - Putamen (basal ganglia)")]="Putamen"
brainregion[which(brainregion=="Brain - Hippocampus")]="Hippocampus"

for (splicetype in splicetypelist){
  #read in the DS analysis result
  inputpath=paste(rootinput,splicetype,counttype,sep="/")
  label=strsplit(rootinput,split="/")[[1]][12]
  pvfilename=paste("pvmatrix_",label,".txt",sep="")
  setwd(inputpath)
  pvmatrix=read.table(pvfilename,sep="\t",header=T)
  
  #get the events information from the row name
  eventsinfo=matrix(NA,dim(pvmatrix)[1],10)
  for (i in 1:dim(pvmatrix)[1]){
    temp=strsplit(rownames(pvmatrix)[i],split="\\|")[[1]]
    eventsinfo[i,]=temp[2:11]
  }
  
  #calculate FDR
  temp=as.numeric(as.matrix(pvmatrix))
  temp[which(is.na(temp))]=1           #assign 1 to missing p values (not converged)
  fdrmatrix=matrix(temp,dim(pvmatrix)[1],dim(pvmatrix)[2])
  rownames(fdrmatrix)=rownames(pvmatrix)
  colnames(fdrmatrix)=colnames(pvmatrix)
  for (i in 1:dim(fdrmatrix)[2]){
    fdrmatrix[,i]=p.adjust(fdrmatrix[,i],"BH")
  }
  
  #get events to plot
  fdrDS.BR=which(fdrmatrix[,"pBR"]<fdrcutoff)    #events pass FDR cutoff
  fdrDS.AGE=which(fdrmatrix[,"pAGE"]<fdrcutoff) 
  fdrDS.GENDER=which(fdrmatrix[,"pGENDER"]<fdrcutoff) 
  
  setwd(outputpath)
  if (length(fdrDS.BR)>0){
    outputtable=cbind(eventsinfo[fdrDS.BR,],pvmatrix[fdrDS.BR,"pBR"],fdrmatrix[fdrDS.BR,"pBR"])
    colnames(outputtable)=c("Exon ID","Gene Symbol","Chr","Strand",
                            "Exon Start","Exon End","Upstream Exon Start","Upstream Exon End","Downstream Exon Start","Downstream Exon End",
                            "P-value","FDR")
    write.table(outputtable,paste(splicetype,"_BR_dependent_splicing_event_summary.txt",sep=""),sep="\t")
  }
  if (length(fdrDS.AGE)>0){
    outputtable=cbind(eventsinfo[fdrDS.AGE,],pvmatrix[fdrDS.AGE,"pAGE"],fdrmatrix[fdrDS.AGE,"pAGE"])
    colnames(outputtable)=c("Exon ID","Gene Symbol","Chr","Strand",
                            "Exon Start","Exon End","Upstream Exon Start","Upstream Exon End","Downstream Exon Start","Downstream Exon End",
                            "P-value","FDR")
    write.table(outputtable,paste(splicetype,"_AGE_dependent_splicing_event_summary.txt",sep=""),sep="\t")
  }
  if (length(fdrDS.GENDER)>0){
    outputtable=cbind(eventsinfo[fdrDS.GENDER,],pvmatrix[fdrDS.GENDER,"pGENDER"],fdrmatrix[fdrDS.GENDER,"pGENDER"])
    colnames(outputtable)=c("Exon ID","Gene Symbol","Chr","Strand",
                            "Exon Start","Exon End","Upstream Exon Start","Upstream Exon End","Downstream Exon Start","Downstream Exon End",
                            "P-value","FDR")
    write.table(outputtable,paste(splicetype,"_GENDER_dependent_splicing_event_summary.txt",sep=""),sep="\t")
  }
}
