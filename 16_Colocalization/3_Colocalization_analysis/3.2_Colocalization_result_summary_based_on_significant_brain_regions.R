###############################
#read in the full result first#
###############################
colocpath="/path/to/14_Colocalization/3_Colocalization_analysis"

setwd(colocpath)
result=read.table("colocalization_result_summary.txt",sep="\t")

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

inputpath="/path/to/14_Colocalization/0_search_and_download_summary_statistics"
exoninfopath="/path/to/input_splicing/logit/JC/SE"

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
  
  rootinput="/path/to/14_Colocalization/2_Colocalization_input/result"   
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


###################################################################
#select the result for sQTL-GWAS pairs in sQTL significant regions#
###################################################################
outputpath="/path/to/coloc/analysis/result"
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



