#Purpose of this code:
#for the sQTL exon we want to do colocalization test on, we get the input (combine sQTL and GWAS information) for colocalization analysis
#So this code will submit each sQTL exon-brain region-gwas study combination as a job, for each job, we get the input

job <- as.numeric(Sys.getenv("SGE_TASK_ID"))
print(paste("job ID:",job))

inputpath="/path/to/14_Colocalization/0_search_and_download_summary_statistics"
rootoutputpath="/path/to/14_Colocalization/3_Colocalization_analysis/result"
exoninfopath="/path/to/input_splicing/logit/JC/SE"
GWASstudyinfopath="/path/to/14_Colocalization/0_search_and_download_summary_statistics"

#read in the information of 36 GWAS studies with summary statistics
setwd(GWASstudyinfopath)
GWASinfo=read.table("Summary_of_downloadable_GWAS_studies.txt",sep="\t",header=T)
rownames(GWASinfo)=GWASinfo[,"rownamelist"]

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

currentcombination=strsplit(combinationlist[job],split="~")[[1]]

testedexon=currentcombination[1]
brainregion=currentcombination[2]
GWASstats=currentcombination[3]

library("coloc")

rootinput="/path/to/14_Colocalization/2_Colocalization_input/result"   ###change###
inputpath=paste(rootinput,gsub("\\|",",",testedexon),GWASstats,brainregion,sep="/")

###################
#read in the input#
###################
setwd(inputpath)
#get the files in the input folder
files=system("ls",intern=T)

if ("no_harmonized_summary_statistics.txt" %in% files){        #if there is no harmonized summary statistics
  print("no_harmonized_summary_statistics")
}

if ("colocalization_input.txt" %in% files){                    #if we have harmonized summary statistics
  colocinput_all=read.table("colocalization_input.txt",sep="\t",header=T)

  #############################
  #run colocalization analysis#
  #############################
  output=matrix(NA,1,6)
  rownames(output)=paste(gsub("\\|",",",testedexon),GWASstats,brainregion,sep="~")
  
  if ("GWAS.n" %in% colnames(colocinput_all)){     #if there is sample size information in the GWAS summary statistics (we need this for colocalization analysis)

    #############################################################################################################################################
    #1. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (original beta from GWAS summary statistics) + MAF from sQTL#
    #############################################################################################################################################
    row2remove=union(union(union(which(is.na(colocinput_all$GWAS.beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(is.na(colocinput_all[,"GWAS.n"]))),which(colocinput_all$GWAS.standard_error==0))     #remove rows with missing value (we cannot do the calculation without these values)
    if (length(row2remove)>0){
      colocinput=colocinput_all[-row2remove,]
    }else{
      colocinput=colocinput_all
    }
    if (dim(colocinput)[1]>0){     #if we have anything left after removing missing values
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(beta=colocinput$GWAS.beta, varbeta=(colocinput$GWAS.standard_error)^2, N=max(colocinput[,"GWAS.n"]),type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,1]=my.res$summary[6]
    }
    
    #######################################################################################################################################
    #2. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (harmonized beta from GWAS statistics) + MAF from sQTL#
    #######################################################################################################################################
    row2remove=union(union(union(which(is.na(colocinput_all$GWAS.hm_beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(is.na(colocinput_all[,"GWAS.n"]))),which(colocinput_all$GWAS.standard_error==0))
    if (length(row2remove)>0){
      colocinput=colocinput_all[-row2remove,]
    }else{
      colocinput=colocinput_all
    }
    if (dim(colocinput)[1]>0){
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(beta=colocinput$GWAS.hm_beta, varbeta=(colocinput$GWAS.standard_error)^2, N=max(colocinput[,"GWAS.n"]),type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,2]=my.res$summary[6]
    }
    
    #########################################################################################################################
    #3. using p value and MAF for sQTL dataset + p value and MAF for GWAS dataset (original GWAS statistics) + MAF from sQTL#
    #########################################################################################################################
    row2remove=which(is.na(colocinput_all[,"GWAS.n"]))
    if (length(row2remove)>0){
      colocinput=colocinput_all[-row2remove,]
    }else{
      colocinput=colocinput_all
    }
    if (dim(colocinput)[1]>0){
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm,N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(pvalues=colocinput$GWAS.p_value,N=max(colocinput[,"GWAS.n"]),type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,3]=my.res$summary[6]
    }
    
    if (sum(is.na(colocinput_all[,"GWAS.hm_effect_allele_frequency"]))<dim(colocinput_all)[1]){         #if the GWAS summary statistics have at least 1 non-NA MAF information
      #############################################################################################################################################
      #4. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (original beta from GWAS summary statistics) + MAF from GWAS#
      #############################################################################################################################################
      row2remove=union(union(union(which(is.na(colocinput_all$GWAS.beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(is.na(colocinput_all[,"GWAS.n"]))),which(colocinput_all$GWAS.standard_error==0))
      if (length(row2remove)>0){
        colocinput=colocinput_all[-row2remove,]
      }else{
        colocinput=colocinput_all
      }
      if (dim(colocinput)[1]>0){
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(beta=colocinput$GWAS.beta, varbeta=(colocinput$GWAS.standard_error)^2, N=max(colocinput[,"GWAS.n"]),type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,4]=my.res$summary[6]
      }
      
      #######################################################################################################################################
      #5. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (harmonized beta from GWAS statistics) + MAF from GWAS#
      #######################################################################################################################################
      row2remove=union(union(union(which(is.na(colocinput_all$GWAS.hm_beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(is.na(colocinput_all[,"GWAS.n"]))),which(colocinput_all$GWAS.standard_error==0))
      if (length(row2remove)>0){
        colocinput=colocinput_all[-row2remove,]
      }else{
        colocinput=colocinput_all
      }
      if (dim(colocinput)[1]>0){
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(beta=colocinput$GWAS.hm_beta, varbeta=(colocinput$GWAS.standard_error)^2, N=max(colocinput[,"GWAS.n"]),type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,5]=my.res$summary[6]
      }
      
      #########################################################################################################################
      #6. using p value and MAF for sQTL dataset + p value and MAF for GWAS dataset (original GWAS statistics) + MAF from GWAS#
      #########################################################################################################################
      row2remove=which(is.na(colocinput_all[,"GWAS.n"]))
      if (length(row2remove)>0){
        colocinput=colocinput_all[-row2remove,]
      }else{
        colocinput=colocinput_all
      }
      if (dim(colocinput)[1]>0){
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm,N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(pvalues=colocinput$GWAS.p_value,N=max(colocinput[,"GWAS.n"]),type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,6]=my.res$summary[6]
      }
    }
    
    setwd(rootoutputpath)
    write.table(output,paste(gsub("\\|",",",testedexon),"~",GWASstats,"~",brainregion,".txt",sep=""),sep="\t",col.names=F)
  }else{         #if there is no sample size information in the summary statistics, we used the one in the GWAS catalog
    currentGWAS.SS=as.numeric(GWASinfo[GWASstats,"samplesize"])      #sample size for the current GWAS study
    print(paste(GWASstats,currentGWAS.SS,sep=" : "))
    #############################################################################################################################################
    #1. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (original beta from GWAS summary statistics) + MAF from sQTL#
    #############################################################################################################################################
    row2remove=union(union(which(is.na(colocinput_all$GWAS.beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(colocinput_all$GWAS.standard_error==0))     #remove rows with missing value (we cannot do the calculation without these values)
    if (length(row2remove)>0){
      colocinput=colocinput_all[-row2remove,]
    }else{
      colocinput=colocinput_all
    }
    if (dim(colocinput)[1]>0){
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(beta=colocinput$GWAS.beta, varbeta=(colocinput$GWAS.standard_error)^2, N=currentGWAS.SS,type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,1]=my.res$summary[6]
    }
    
    #######################################################################################################################################
    #2. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (harmonized beta from GWAS statistics) + MAF from sQTL#
    #######################################################################################################################################
    row2remove=union(union(which(is.na(colocinput_all$GWAS.hm_beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(colocinput_all$GWAS.standard_error==0))
    if (length(row2remove)>0){
      colocinput=colocinput_all[-row2remove,]
    }else{
      colocinput=colocinput_all
    }
    if (dim(colocinput)[1]>0){
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(beta=colocinput$GWAS.hm_beta, varbeta=(colocinput$GWAS.standard_error)^2, N=currentGWAS.SS,type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,2]=my.res$summary[6]
    }
    
    #########################################################################################################################
    #3. using p value and MAF for sQTL dataset + p value and MAF for GWAS dataset (original GWAS statistics) + MAF from sQTL#
    #########################################################################################################################
    if (length(which(is.na(colocinput_all$GWAS.p_value)))<dim(colocinput_all)[1]){      #if we have at least one non-missing p value from GWAS summary statistics (there is one exon that all the GWAS p values are missing)
      colocinput=colocinput_all
      my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm,N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                          dataset2=list(pvalues=colocinput$GWAS.p_value,N=currentGWAS.SS,type="quant"),
                          MAF=colocinput$sQTL.MAF)
      output[1,3]=my.res$summary[6]
    }

    if (sum(is.na(colocinput_all[,"GWAS.hm_effect_allele_frequency"]))<dim(colocinput_all)[1]){         #if the GWAS summary statistics have at least 1 non-NA MAF information
      #############################################################################################################################################
      #4. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (original beta from GWAS summary statistics) + MAF from GWAS#
      #############################################################################################################################################
      row2remove=union(union(which(is.na(colocinput_all$GWAS.beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(colocinput_all$GWAS.standard_error==0))
      if (length(row2remove)>0){
        colocinput=colocinput_all[-row2remove,]
      }else{
        colocinput=colocinput_all
      }
      if (dim(colocinput)[1]>0){
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(beta=colocinput$GWAS.beta, varbeta=(colocinput$GWAS.standard_error)^2, N=currentGWAS.SS,type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,4]=my.res$summary[6]
      }
      
      #######################################################################################################################################
      #5. p value and MAF for sQTL dataset + beta and standard error for GWAS dataset (harmonized beta from GWAS statistics) + MAF from GWAS#
      #######################################################################################################################################
      row2remove=union(union(which(is.na(colocinput_all$GWAS.hm_beta)),which(is.na(colocinput_all$GWAS.standard_error))),which(colocinput_all$GWAS.standard_error==0))
      if (length(row2remove)>0){
        colocinput=colocinput_all[-row2remove,]
      }else{
        colocinput=colocinput_all
      }
      if (dim(colocinput)[1]>0){
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm, N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(beta=colocinput$GWAS.hm_beta, varbeta=(colocinput$GWAS.standard_error)^2, N=currentGWAS.SS,type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,5]=my.res$summary[6]
      }
      
      #########################################################################################################################
      #6. using p value and MAF for sQTL dataset + p value and MAF for GWAS dataset (original GWAS statistics) + MAF from GWAS#
      #########################################################################################################################
      if (length(which(is.na(colocinput_all$GWAS.p_value)))<dim(colocinput_all)[1]){      #if we have at least one non-missing p value from GWAS summary statistics (there is one exon that all the GWAS p values are missing)
        colocinput=colocinput_all
        my.res <- coloc.abf(dataset1=list(pvalues=colocinput$sQTL.pvals.lm,N=max(colocinput[,"sQTL.Sample.size"]),type="quant"),
                            dataset2=list(pvalues=colocinput$GWAS.p_value,N=currentGWAS.SS,type="quant"),
                            MAF=colocinput$GWAS.hm_effect_allele_frequency)
        output[1,6]=my.res$summary[6]
      }
    }
    
    setwd(rootoutputpath)
    write.table(output,paste(gsub("\\|",",",testedexon),"~",GWASstats,"~",brainregion,".txt",sep=""),sep="\t",col.names=F)
  }
}


