#the purpose of this code is to check the total number of colocalization analysis, i.e., total number of expected output files

library("coloc")

inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis/result"
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
      
      singlecombination=paste(exonfullID,gwasstudy,currentBR,sep="~")
      
      combinationlist=c(combinationlist,singlecombination)
    }
  }
}

combinationlist=unique(combinationlist)        #there are duplications  


totalrun=0       #total number of colocalizaiton calculation
outputfilelist=c()
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
    
    filename=paste(gsub("\\|",",",testedexon),"~",GWASstats,"~",brainregion,".txt",sep="")
    outputfilelist=c(outputfilelist,filename)
    
    if ("GWAS.n" %in% colnames(colocinput_all)){     #if there is sample size information in the GWAS summary statistics (we need this for colocalization analysis)
      totalrun=totalrun+1
    }else{
       totalrun=totalrun+1                           #now we can run without sample size information in the summary statistics
     }
  }
}


totalrun     #number of runs when we require sample size information in summary statistics
#[1] 130
totalrun     #number of runs when we don't require sample size information in summary statistics
#[1] 559
length(unique(outputfilelist))
#559



#######################################################################
#check missing files and their job ID when the calculation is finished#
#######################################################################
setwd("/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/3_Colocalization_analysis/result")
command="ls *.txt"
outputlist=system(command,intern=T)

#missing files:
missingfile=setdiff(outputfilelist,outputlist)

#missing file job ID:
missingID=c()
for (i in 1:length(missingfile)){
  missingID=c(missingID,which(combinationlist %in% strsplit(gsub(",","\\|",missingfile[i]),split=".txt")[[1]][1]))
}


