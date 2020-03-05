#Purpose of this code:
#for the sQTL exon we want to do colocalization test on, we get the input (combine sQTL and GWAS information) for colocalization analysis
#So this code will submit each sQTL exon-brain region-gwas study combination as a job, for each job, we get the input

###############################################################get the exon-brain region-gwas study combination####################################################################################
job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(paste("job ID:",job))


inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/2_Colocalization_input/result"
exoninfopath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE"


setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#exon full ID - short ID conversion
setwd(exoninfopath)
exonfullinfo=read.table("exon_info.fromGTF.SE.txt",sep="\t",header=T)
rownames(exonfullinfo)=exonfullinfo[,"ID"]

#brain region name conversion
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
formalbrainregionlist=c("Amygdala",
                        "Anterior cingulate cortex BA24",
                        "Caudate basal ganglia",
                        "Cerebellar Hemisphere",
                        "Cerebellum",
                        "Cortex",
                        "Frontal Cortex BA9",
                        "Hippocampus",
                        "Hypothalamus",
                        "Nucleus accumbens basal ganglia",
                        "Putamen basal ganglia",
                        "Spinal cord cervical c-1",
                        "Substantia nigra")


#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)


#get all the exon-brain region-GWAS study combinations
combinationlist=c()

for (i in 1:dim(gwasftp)[1]){                      #for each exon
  exonshortID=as.character(gwasftp[i,"exon_short_ID"])
  
  #get the full exon ID
  exonfullID=paste(c(strsplit(exonshortID,split="_")[[1]][2],as.matrix(exonfullinfo[exonshortID,])[-1]),collapse="|")
  
  #get all the significant brain regions
  #sig.region=strsplit(uniqueallsigregion[i],split=", ")[[1]]
  sig.region=formalbrainregionlist         #we carry out the calculation on all brain regions, not just the significant brain regions
  
  #get all the GWAS studies
  firstauthorlist=strsplit(as.character(gwasftp[i,"First.Author"]),split=",")[[1]][-1]
  pubmedidlist=strsplit(as.character(gwasftp[i,"PubMed.ID"]),split=",")[[1]][-1]
  studyaccesslist=strsplit(as.character(gwasftp[i,"Study.accession"]),split=",")[[1]][-1]
  ftplinklist=strsplit(as.character(gwasftp[i,"FTP.link"]),split=",")[[1]][-1]
  #remove information from GWAS studies with no FTP link (downloadable=1 means at least one GWAS study has downloadable summary statistics, not all of them)
  firstauthorlist=firstauthorlist[which(ftplinklist!="no_link")]
  pubmedidlist=pubmedidlist[which(ftplinklist!="no_link")]
  studyaccesslist=studyaccesslist[which(ftplinklist!="no_link")]
  
  #get the significant brain regions and generate folder for each brain region
  for (br in 1:length(sig.region)){
    currentBR=brainregionlist[which(formalbrainregionlist %in% sig.region[br])]
    
    #get the corresponding GWAS study
    for (j in 1:length(firstauthorlist)){
      gwasstudy=paste(firstauthorlist[j],pubmedidlist[j],studyaccesslist[j],sep="_")
      
      singlecombination=paste(exonfullID,currentBR,gwasstudy,sep="~")
      
      combinationlist=c(combinationlist,singlecombination)
    }
  }
}

combinationlist=unique(combinationlist)        #there are duplications  
      
 
###############################################################start the calculation for each combination####################################################################################
currentcombination=strsplit(combinationlist[job],split="~")[[1]]

testedexon=currentcombination[1]
brainregion=currentcombination[2]
GWASstats=currentcombination[3]

exonbrgwaspath=paste(rootoutputpath,gsub("\\|",",",testedexon),GWASstats,brainregion,sep="/")
command=paste("mkdir -p",exonbrgwaspath)
system(command)


library("data.table")

rootinput="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/1_MAF_calculation"             #input from sQTL analysis
inputpath=paste(rootinput,"/result","/",gsub("\\|",",",testedexon),"/",brainregion,sep="")

rootgwassspath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"

rootoutput="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/2_Colocalization_input/result"
outputpath=paste(rootoutput,"/",gsub("\\|",",",testedexon),"/",GWASstats,"/",brainregion,sep="")

#############################
#read in the MAF information#
#############################
setwd(inputpath)
MAF.SS=read.table("Colocalization_input.txt",sep="\t",header=T)

#remove rows with no rsID
if (length(which(MAF.SS[,"rsID"]=="."))>0){           #if we have SNPs with no rsID
  MAF.SS=MAF.SS[-which(MAF.SS[,"rsID"]=="."),]
}
rownames(MAF.SS)=MAF.SS[,"rsID"]

#####################################
#read in the GWAS summary statistics#
#####################################
gwaspath=paste(rootgwassspath,GWASstats,"harmonised",sep="/")
#check if there is harmonized data (this pipeline only applies to GWAS studies with harmonized summary statistics)
if (dir.exists(gwaspath)){
  setwd(gwaspath)
  
  #get the file names
  filelist=system("ls",intern=T)
  
  #we use the harmonised_file
  GWAS <-fread(filelist[grep("h.tsv",filelist)], header=T, sep="\t") 
  GWAS=as.data.frame(GWAS)
  
  #############################################################
  #get the overlap between the two and generate the input file#
  #############################################################
  sQTLrsID=as.character(MAF.SS[,"rsID"])
  
  GWASrsID1=as.character(GWAS[,"hm_rsid"])
  GWASrsID2=as.character(GWAS[,"variant_id"])
  GWASrsID=GWASrsID1
  missingrsIDpos=which(is.na(GWASrsID))
  GWASrsID[missingrsIDpos]=GWASrsID2[missingrsIDpos]      #for missing rsID in GWASrsID1, we replace with information in GWASrsID2
    
  overlaprsID=intersect(sQTLrsID,GWASrsID)
  print(paste("number of overlapped SNPs:",length(overlaprsID)))
  
  if (length(overlaprsID)>0){
    output=matrix(NA,length(overlaprsID),dim(GWAS)[2]+dim(MAF.SS)[2])
    colnames(output)=c(paste("sQTL.",colnames(MAF.SS),sep=""),paste("GWAS.",colnames(GWAS),sep=""))
    rownames(output)=overlaprsID
    
    for (i in 1:length(overlaprsID)){
      print(i)
      infofromsQTL=MAF.SS[overlaprsID[i],]
      infofromGWAS=GWAS[which(GWAS[,"hm_rsid"] %in% overlaprsID[i])[1],]      #there may be multiple hits (some of them are duplicated, some are not, we just choose the first one)
      output[i,1:dim(MAF.SS)[2]]=as.character(as.matrix(infofromsQTL))
      output[i,(dim(MAF.SS)[2]+1):dim(output)[2]]=as.character(as.matrix(infofromGWAS))
    }
    
    setwd(outputpath)
    write.table(output,"colocalization_input.txt",sep="\t")
  }else{
    setwd(outputpath)
    info=paste(brainregion,testedexon,GWASstats,sep="\t")
    write.table(info,"no_overlapped_rsID.txt",sep="\t",row.names=F,col.names=F,quote=F)
  }
}else{
  setwd(outputpath)
  info=paste(brainregion,testedexon,GWASstats,sep="\t")
  write.table(info,"no_harmonized_summary_statistics.txt",sep="\t",row.names=F,col.names=F,quote=F)
}


















