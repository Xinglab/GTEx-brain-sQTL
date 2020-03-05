#Purpose of this code:
#For all the 36 GWAS studies with downloadable summary statistics, we collect their information from GWAS catalog
#we also annotate the 36 GWAS studies to show if they have harmonized summary statistics and if they have sample size information in their summary statistics

outputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
SSfolder="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/GWAS_summary_statistics"

#read in the GWAS catalog table#
GWASdb="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/GWAS_databases/gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
GWAS=read.table(GWASdb,sep="\t",header=T,fill=TRUE,check.names=F,quote="",comment.char="")

#read in the list of downloadable summary statistics
inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/14_Colocalization/0_search_and_download_summary_statistics"
setwd(inputpath)
gwas.ss=read.table("sQTL_GWAS_summary_statistics_downloadable_SE_logit_JC_pvalue.txt",sep="\t",header=T)


#get the GWAS information of the downloadable GWAS studies
downloadableGWAS=matrix(NA,0,dim(GWAS)[2])
colnames(downloadableGWAS)=colnames(GWAS)
#add additional information
additionalinfo=matrix(NA,0,2)
colnames(additionalinfo)=c("have_harmonized_summary_statistics?","have_sample_size_for_colocalization_analysis?")

rownamelist=c()

for (i in 1:dim(gwas.ss)[1]){
  authorlist=strsplit(as.character(gwas.ss[i,"First.Author"]),split=",")[[1]][-1]
  pubmedidlist=strsplit(as.character(gwas.ss[i,"PubMed.ID"]),split=",")[[1]][-1]
  studyaccesslist=strsplit(as.character(gwas.ss[i,"Study.accession"]),split=",")[[1]][-1]
  ftplinklist=strsplit(as.character(gwas.ss[i,"FTP.link"]),split=",")[[1]][-1]
  
  for (j in 1:length(ftplinklist)){
    if (ftplinklist[j]=="no_link"){      #for sQTLs with multiple GWAS studies, not all of them have summary statistics to download
      donothing="donothing"
    }else{
      author=authorlist[j]
      pubmedid=pubmedidlist[j]
      studyaccess=studyaccesslist[j]
      ftplink=ftplinklist[j]
      
      #get the GWAS information of this study
      rownum=intersect(which(GWAS[,"PUBMEDID"] %in% pubmedid),which(GWAS[,"STUDY ACCESSION"] %in% studyaccess))
      col2select=c("DATE ADDED TO CATALOG",
                   "PUBMEDID",
                   "FIRST AUTHOR",
                   "DATE",
                   "JOURNAL",
                   "LINK",
                   "STUDY",
                   "DISEASE/TRAIT",
                   "INITIAL SAMPLE SIZE",        
                   "REPLICATION SAMPLE SIZE",
                   "MAPPED_TRAIT",
                   "MAPPED_TRAIT_URI",          
                   "STUDY ACCESSION")      #we only care about information about the GWAS study, however, GWAS catalog contains also the SNPs reported by each study, we don't need this information
      currentGWASstudy=unique(GWAS[rownum,col2select])
      
      #add the information to the result
      downloadableGWAS=rbind(downloadableGWAS,currentGWASstudy)
      
      #check if this study has harnomized summary statistics or not & have sample size
      label=paste(author,pubmedid,studyaccess,sep="_")
      rownamelist=c(rownamelist,label)
      setwd(paste(SSfolder,label,sep="/"))
      
      if (dir.exists("harmonised")){        #if there is harmonized data
          harmonized="Yes"
          
          #check if there is sample size information
          setwd(paste(SSfolder,label,"harmonised",sep="/"))
          #get the file names
          filelist=system("ls",intern=T)
          
          #we use the harmonised_file (get the header information, we don't need the whole file)
          command=paste("head -1",filelist[grep("h.tsv",filelist)])
          SS.header <- system(command,intern=T)
          
          if ("n" %in% strsplit(SS.header,split="\t")[[1]]){      #if the summary statistics have sample size information
            samplesizeinfo="Yes"
          }else{
            samplesizeinfo="No"
          }
      }else{
        harmonized="No"
        samplesizeinfo="No"
      }
      
      #add the information to the result
      temp=matrix(NA,1,2)
      temp[1,1]=harmonized
      temp[1,2]=samplesizeinfo
      additionalinfo=rbind(additionalinfo,temp)
    }
  }
}

output=cbind(downloadableGWAS,additionalinfo,rownamelist)
output=unique(output)

#get the sample size information of each GWAS study
samplesize=rep(NA,dim(output)[1])
for (i in 1:dim(output)[1]){
  rawinfo=paste(as.character(output[i,"INITIAL SAMPLE SIZE"]),as.character(output[i,"REPLICATION SAMPLE SIZE"]),sep=" ")     
  #we use the information in the two columns to calculate sample size
  #reason: for GWAS studies with summary statistics, I compared the sample size there and the sample size here. 
  #If we use the information in just one column, the sample size calculated by information in the GWAS catalog can be smaller than the sample size reported in the summary statistics, which should not happen
  #so we need to add up all the samples in the two columns
  #get the numbers from this string
  #1. replace "," with ""
  rawinfo=gsub(",", "",  rawinfo)
  #2. replace non-numeric values to "-"
  rawinfo=gsub("[^0-9]", " ",  rawinfo)     #if we have decimal point, we need to use gsub("[^0-9.]", " ",  rawinfo), but here we have all integers so we don't need it  
  temp=unique(strsplit(rawinfo,split=" ")[[1]])
  samplesize[i]=sum(as.numeric(temp[-which(temp=="")]))
}

samplesize[19]=samplesize[19]-508

output=cbind(output,samplesize)

setwd(outputpath)
write.table(output,"Summary_of_downloadable_GWAS_studies.txt",sep="\t")

