#Purpose of this code:
#Check the total number of GWAS studies associated with GWAS associated sQTLs

#########################
#read in the sQTL result#
#########################
splicetype="SE"
PSItype="logit"
counttype="JC"
type="pvalue"
rootinput="/path/to/summary"

setwd(paste(rootinput,PSItype,counttype,splicetype,sep="/"))
sqtl=read.table(paste("sQTL_summary_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type,"_sorted.txt",sep=""),sep="\t",header=T,quote="")

################################
#read in the GWAS catalog table#
################################
GWASdb="/path/to/GWAS_catalog/files/gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
GWAS=read.table(GWASdb,sep="\t",header=T,fill=TRUE,check.names=F,quote="",comment.char="")

disease_key_word=toupper(c("Alzheimer","Amyotrophic lateral sclerosis","Parkinson","frontotemporal dementia","Huntington",
                           "epilepsy","autism","schizophrenia","bipolar","depression",
                           "attention deficit hyperactivity disorder","glio","multiple sclerosis","narcolepsy","stroke"))       #change all the names to upper case for comparison


####################################################
#1. get sQTLs with GWAS association
sqtl_select=subset(sqtl,!is.na(sqtl[,"disease.ontology"]))

#2. get their GWAS study information for each sQTL 
GWAS_SS=c()

for (i in 1:dim(sqtl_select)[1]){        #for each sQTL

  #for each sQTL event, get the GWAS SNP rsID
  rsID=as.character(sqtl_select[i,"GWAS.rsID"])
  
  #get the GWAS study pubmed ID(s) and study accession number(s) of the GWAS snp
  info=GWAS[which(GWAS[,"SNPS"] %in% rsID),]

  for (j in 1:dim(info)[1]){       #there may be multiple GWAS studies related to the GWAS SNP
    temppubmedID=info[j,"PUBMEDID"]
    tempstudyaccess=as.character(info[j,"STUDY ACCESSION"])
    tempfirstauthor=gsub(" ", "", as.character(info[j,"FIRST AUTHOR"]), fixed = TRUE) 
    
    uniqueID=paste(tempfirstauthor,temppubmedID,tempstudyaccess,sep="_")
    
    GWAS_SS=c(GWAS_SS,uniqueID)
  }
}

length(unique(GWAS_SS))
#278



