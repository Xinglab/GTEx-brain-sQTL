#the purpose of this code is to summarize all the colocalization result
#The difference between this code and the result in /u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis is that:
#For the colocalization calculation, we calculate that for each sQTL-GWAS pair in all 13 studies. 
#But here, we want to summarize the result based on those pairs in sQTL significant regions instead of all 13 regions.

###############################
#read in the full result first#
###############################
colocpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis"

setwd(colocpath)
result=read.table("colocalization_result_summary.txt",sep="\t")

#rownames(result)=result[,1]
#result=result[,-1]

colnames(result)=c("sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF","sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF","sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF",
                   "sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF","sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF","sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF")


  
##############################################################
#get the list and number of colocalization analysis performed#
##############################################################
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

inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
exoninfopath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE"

setwd(inputpath)
GWAS_SS=read.table("sQTL_GWAS_summary_statistics_SE_logit_JC_pvalue.txt",sep="\t",header=T)

#get the subset with downloadable summary statistics
gwasftp=subset(GWAS_SS,GWAS_SS[,"downloadable"]==1)

#get the list of unique exons
exonlist=as.character(unique(gwasftp[,"exon_short_ID"]))

#get the significant brain regions of each exon (for the same exon, the SNP could be different so the significant brain regions could also be different, here we take the union for each exon)
uniqueallsigregion=rep(NA,length(exonlist))
for (i in 1:length(exonlist)){
  exonshortID=exonlist[i]
  allsigregion=as.character(gwasftp[which(gwasftp[,"exon_short_ID"] %in% exonshortID),"significant.region"])
  if (length(allsigregion)==1){           #if the exon only appears once
    uniqueallsigregion[i]=allsigregion
  }else{                                  #if the exon appears multiple times (either because multiple GWAS studies or different sQTL SNPs)
    uniqueallsigregion[i]=paste(unique(strsplit(paste(allsigregion,collapse=", "),split=", ")[[1]]),collapse=", ")
  }
}

#exon full ID - short ID conversion
setwd(exoninfopath)
exonfullinfo=read.table("exon_info.fromGTF.SE.txt",sep="\t",header=T)
rownames(exonfullinfo)=exonfullinfo[,"ID"]

#get all the exon-brain region-GWAS study combinations (this contains all the GWAS studies with downloadable summary statistics, but some of them may not have harmoznied summary statistics)
combinationlist=c()

for (i in 1:dim(gwasftp)[1]){                      #for each exon
  exonshortID=as.character(gwasftp[i,"exon_short_ID"])
  
  #get the full exon ID
  exonfullID=paste(c(strsplit(exonshortID,split="_")[[1]][2],as.matrix(exonfullinfo[exonshortID,])[-1]),collapse="|")
  
  #get all the significant brain regions of the current sQTL event
  sig.region=strsplit(uniqueallsigregion[which(exonlist %in% exonshortID)],split=", ")[[1]]
  
  #get all the GWAS studies
  firstauthorlist=strsplit(as.character(gwasftp[i,"First.Author"]),split=",")[[1]][-1]
  pubmedidlist=strsplit(as.character(gwasftp[i,"PubMed.ID"]),split=",")[[1]][-1]
  studyaccesslist=strsplit(as.character(gwasftp[i,"Study.accession"]),split=",")[[1]][-1]
  ftplinklist=strsplit(as.character(gwasftp[i,"FTP.link"]),split=",")[[1]][-1]
  #remove information from GWAS studies with no FTP link (downloadable=1 means at least one GWAS study has downloadable summary statistics, not all of them)
  firstauthorlist=firstauthorlist[which(ftplinklist!="no_link")]
  pubmedidlist=pubmedidlist[which(ftplinklist!="no_link")]
  studyaccesslist=studyaccesslist[which(ftplinklist!="no_link")]
  
  #get the significant brain regions
  for (br in 1:length(sig.region)){
    currentBR=brainregionlist[which(formalbrainregionlist %in% sig.region[br])]
    
    #get the corresponding GWAS study
    for (j in 1:length(firstauthorlist)){
      gwasstudy=paste(firstauthorlist[j],pubmedidlist[j],studyaccesslist[j],sep="_")
      
      singlecombination=paste(exonfullID,gwasstudy,currentBR,sep="~")
      
      combinationlist=c(combinationlist,singlecombination)
    }
  }
}

combinationlist=unique(combinationlist)        #there are duplications  


totalrun=0       #total number of colocalizaiton calculation
outputfilelist=c()        #list of exon-brain region-GWAS study combinations (we only include those combinations containing GWAS studies with harmonized summary statistics)
for (job in 1:length(combinationlist)){
  currentcombination=strsplit(combinationlist[job],split="~")[[1]]
  
  testedexon=currentcombination[1]
  brainregion=currentcombination[3]
  GWASstats=currentcombination[2]
  
  rootinput="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/2_Colocalization_input/result"   
  inputpath=paste(rootinput,gsub("\\|",",",testedexon),GWASstats,brainregion,sep="/")
  
  ###################
  #read in the input#
  ###################
  setwd(inputpath)
  #get the files in the input folder
  files=system("ls",intern=T)
  
  if ("colocalization_input.txt" %in% files){                    #if we have harmonized summary statistics
    colocinput_all=read.table("colocalization_input.txt",sep="\t",header=T)
    
    filename=paste(gsub("\\|",",",testedexon),"~",GWASstats,"~",brainregion,sep="")
    outputfilelist=c(outputfilelist,filename)
    
    totalrun=totalrun+1  
  }
}

totalrun     #number of runs when we don't require sample size information in summary statistics
#[1] 124
length(unique(outputfilelist))
#124


###################################################################
#select the result for sQTL-GWAS pairs in sQTL significant regions#
###################################################################
outputpath="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/colocalization_analysis_summary_based_on_sQTL_significant_brain_regions/result"
result.sigregion=result[unique(outputfilelist),]

#change the format of the result table
combination.info=matrix(NA,dim(result.sigregion)[1],5)
colnames(combination.info)=c("Exon",	"GWAS study - First author",	"GWAS study - Pubmed ID",	"GWAS study - Study accession",	"Brain region")
for (i in 1:dim(result.sigregion)[1]){
  temp=strsplit(rownames(result.sigregion)[i],split="~")[[1]]
  gwas.temp=strsplit(temp[2],split="_")[[1]]
  combination.info[i,"Exon"]=temp[1]
  combination.info[i,"GWAS study - First author"]=gwas.temp[1]
  combination.info[i,"GWAS study - Pubmed ID"]=gwas.temp[2]
  combination.info[i,"GWAS study - Study accession"]=gwas.temp[3]
  combination.info[i,"Brain region"]=temp[3]
}
result.sigregion=cbind(combination.info,result.sigregion)

setwd(outputpath)
write.table(result.sigregion,"colocalization_result_summary_based_on_sig_region.txt",sep="\t")


########################################
#check the number of significant result#
########################################
inputpath="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/colocalization_analysis_summary_based_on_sQTL_significant_brain_regions/result"
setwd(inputpath)
result=read.table("colocalization_result_summary_based_on_sig_region.txt",sep="\t",header=T)

for (i in 6:11){
  for (j in c(0.75,0.8)){
    print(paste(sum(result[,i]>=j,na.rm=T),"events with PP4 greater or equal to",j,"based on",colnames(result)[i]))
  }
}


##################################################
#not requiring_sample_size_in_summary_statistics:#
##################################################
[1] "56 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"
[1] "46 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.beta.se_sQTL.MAF"

[1] "56 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"
[1] "46 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.harmonized.beta.se_sQTL.MAF"

[1] "77 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"
[1] "63 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.p.MAF_sQTL.MAF"

[1] "17 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF"
[1] "14 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.beta.se_GWAS.MAF"

[1] "17 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF"
[1] "14 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.harmonized.beta.se_GWAS.MAF"

[1] "26 events with PP4 greater or equal to 0.75 based on sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF"
[1] "22 events with PP4 greater or equal to 0.8 based on sQTL.p.MAF_GWAS.original.p.MAF_GWAS.MAF"


#how many sQTL-GWAS pairs
inputpath="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/colocalization_analysis_summary_based_on_sQTL_significant_brain_regions/result"
setwd(inputpath)
result=read.table("colocalization_result_summary_based_on_sig_region.txt",sep="\t",header=T)

length(unique(paste(result[,1],paste(result[,2],result[,3],result[,4],sep="_"),sep="~")))
