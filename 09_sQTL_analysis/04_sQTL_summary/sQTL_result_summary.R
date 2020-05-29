splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
counttype="JC"       #JCEC (junction count + exon body count) or JC (junction count only)
PSItype="logit"      #logit (logit transformed PSI value) or original (original PSI from 0 to 1)
exoninfopath=paste("/path/to/input_splicing",
                   PSItype,counttype,splicetype,sep="/")
codepath="/path/to/sQTL/calculation/code"
rootoutput=paste("/path/to/sQTL_run",
                 PSItype,counttype,splicetype,sep="/")
outputpath=paste("/output/path/summary",PSItype,counttype,splicetype,sep="/")


##################################
#1. exon-sQTL-GWAS result summary#
##################################
type=c("pvalue","permutation")
sQTLLDGWASfolder=c(paste("selected_sQTL_LD_GWAS_",type[1],sep=""),paste("selected_sQTL_LD_GWAS_",type[2],sep=""))

#read in exon information file
setwd(exoninfopath)
exoninfo=read.table(paste("exon_info.fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)
rownames(exoninfo)=exoninfo[,"ID"]
exonshorterID=rownames(exoninfo)

#read in brain region file
setwd(codepath)
brainregion=as.character(as.matrix(read.table("brainregionlist.txt",sep="\t")))

#get genomic position based on SNP information
get_genopos=function(pos){
  temp=strsplit(pos,split="_")[[1]]
  return(paste(temp[1],temp[2],sep="_"))
}


for (j in 1:2){
  
  disease_key_word=toupper(c("Alzheimer","Amyotrophic lateral sclerosis","Parkinson","frontotemporal dementia","Huntington",
                             "epilepsy","autism","schizophrenia","bipolar","depression",
                             "attention deficit hyperactivity disorder","glio","multiple sclerosis","narcolepsy","stroke"))       #change all the names to upper case for comparison
  #summarize how many unique GWAS loci are related to the diseases we are interested in for each brain region
  diseaseGWAS=matrix(0,length(disease_key_word),length(brainregion)) 
  rownames(diseaseGWAS)=disease_key_word
  colnames(diseaseGWAS)=brainregion

  #updated version 
  #(if the sQTL in region A is not the sQTL in region B (i.e., the SNP is the top SNP for an exon in region A but not the top SNP for the same exon in region B)
  #but the p value of this exon-SNP pair is smaller than the p value cutoff of region B, 
  #then the GWAS loci that are in high LD with this sQTL/SNP are not only counted in region A but also counted in region B)
  updated_diseaseGWAS=matrix(0,length(disease_key_word),length(brainregion))    
  rownames(updated_diseaseGWAS)=disease_key_word
  colnames(updated_diseaseGWAS)=brainregion
    
  #summarize how many unique exons are related to the diseases we are interested in for each brain region
  diseaseexon=matrix(0,length(disease_key_word),length(brainregion)) 
  rownames(diseaseexon)=disease_key_word
  colnames(diseaseexon)=brainregion
  
  #updated version 
  #(if the sQTL in region A is not the sQTL in region B (i.e., the SNP is the top SNP for an exon in region A but not the top SNP for the same exon in region B)
  #but the p value of this exon-SNP pair is smaller than the p value cutoff of region B, 
  #then the exon of this sQTL/SNP are not only counted in region A but also counted in region B)
  updated_diseaseexon=matrix(0,length(disease_key_word),length(brainregion)) 
  rownames(updated_diseaseexon)=disease_key_word
  colnames(updated_diseaseexon)=brainregion
  
  totalsqtlnum=rep(0,length(brainregion))       #Total number of significant sQTLs identified in each brain region regardless of GWAS signal (sQTL here means sQTL SNPs, not the exon) and rsID (some of the sQTLs here don't have corresponding rsID)
  totalsqtlexonnum=rep(0,length(brainregion))   #Total number of exons with significant sQTL SNPs regardless of GWAS signal and rsID
  sqtlnum=rep(0,length(brainregion))       #Total number of significant sQTLs identified in each brain region AND correlated with at least one GWAS loci (sQTL here means sQTL SNPs, not the exon)
  GWASnum=rep(0,length(brainregion))       #Total number of unique GWAS loci that are in high LD with the significant sQTLs in each brain region 
  diseaseGWASperregion=rep(0,length(brainregion))   #Total number of unique GWAS loci that are related to at least one disease we are interested in in each brain region
  diseaseexonperregion=rep(0,length(brainregion))   #Total number of unique exons that are related to at least one disease we are interested in in each brain region
  
  numGWASperdisease=list()      #list of unique GWAS loci for each disease across all brain regions
  for (i in disease_key_word){
    numGWASperdisease[i]=i    #we put the name of the disease as a place holder, it needs to be removed when we count the number later
  }
  
  numexonperdisease=list()      #list of unique exons for each disease across all brain regions
  for (i in disease_key_word){
    numexonperdisease[i]=i    #we put the name of the disease as a place holder, it needs to be removed when we count the number later
  }
  
  all_output=matrix(,nrow=0, ncol=20)
  
  for (i in 1:length(brainregion)){
    br=brainregion[i]
    inputpath=paste(rootoutput,"/",br,"/",sQTLLDGWASfolder[j],sep="")
    setwd(inputpath)
    
    LDresult=try(suppressMessages(read.table("sig_SNPS_gwas.ld.LD.result.txt",sep="\t",fill=TRUE,quote=NULL)),silent=TRUE)   #the same sQTL can have multiple high LD GWAS loci which means that there are duplications in sQTL ID
    
    if (!(inherits(LDresult,"try-error"))){       #if the LDresult is not empty (there are significant sQTLs)
      colnames(LDresult)=c("sQTL","highLD_GWAS","disease.trait","disease.ontology","p.value","PUBMEDID","r.square")
      #disease.trait corresponds to the DISEASE/TRAIT column in the GWASdb table
      #disease.ontology corresponds to the MAPPED_TRAIT column in the GWASdb table
      
      #how many sQTL & GWAS loci for each brain region?
      #print(paste(br,dim(LDresult)[1],sep=" : "))
      sqtlnum[i]=length(unique(LDresult[,"sQTL"]))
      GWASnum[i]=length(unique(LDresult[,"highLD_GWAS"]))
      
      #how many unique GWAS loci for each disease in each brain region?
      SNPtrait=toupper(LDresult[,"disease.ontology"])
      for (disease in disease_key_word){
        uniqueGWASID=as.character(unique(LDresult[which(grepl(disease,SNPtrait)),"highLD_GWAS"]))
        diseaseGWAS[disease,br]=length(uniqueGWASID)
        numGWASperdisease[[disease]]=c(numGWASperdisease[[disease]],uniqueGWASID)   #add those unique GWAS traits in the current brain region for each disease
      }
      
      ###########################################################################################
      #connect each exon to its sQTL and then to the GWAS loci that is in high LD with that sQTL#
      ###########################################################################################
      #columns in each brain region
      #exon_fullI_D1  exon_short_ID genomic_position sQTL_rsID  +related_information  GWAS_high_LD_SNP_rsID  +related_information  related_to_disease_or_not high/low/no_significant_in_current_brain_region high/low/no_significant_in_age exon_length 
      
      #exon full ID -> exon shorter ID -> genomic position -> sQTL rsID -> GWAS rsID
      #exon full ID -> exon shorter ID
      exonfullID=rep(NA,length(exonshorterID))
      for (id in 1:length(exonfullID)){
        temp=exoninfo[id,]
        temp[1]=strsplit(as.character(exoninfo[id,"ID"]),split="_")[[1]][2]
        exonfullID[id]=paste(as.matrix(temp),collapse="|")
      }
      output=cbind(exonfullID,exonshorterID)
      rownames(output)=exonshorterID
      colnames(output)=c("exon_full_ID","exon_short_ID")
      
      #exon shorter ID -> genomic position
      inputpath=paste(rootoutput,"/",br,sep="")
      setwd(inputpath)
      sqtlinfo=read.table(paste("selected.sQTL.glmm.Glimmps_each_exon_cis_",br,".",type[j],".txt",sep=""),sep="\t")      #sQTL result (exon-SNP pair)
      totalsqtlnum[i]=length(unique(sqtlinfo[,3]))    #one SNP can be the sQTL of multiple exons
      totalsqtlexonnum[i]=length(unique(sqtlinfo[,1]))
      genoposall=sapply(as.character(sqtlinfo[,3]),get_genopos)
      sqtlinfo=cbind(sqtlinfo,genoposall)
      rownames(sqtlinfo)=sqtlinfo[,1]
      colnames(sqtlinfo)=c("exon_short_ID","chr","SNP_info","relative_pos","sQTL_p_value","sQTL_beta","absolute_distance","relative_distance","genoposall")
      output=subset(output,rownames(output) %in% sqtlinfo[,1])     #filter out exons without significant sQTLs
      output=cbind(output,sqtlinfo[rownames(output),])        #combine the two informaiton
      output=output[,-3]   #remove duplicated columns
      
      #genomic position -> sQTL rsID
      IDmap=read.table(paste("selected.sQTL.glmm.Glimmps_each_exon_cis_",br,".",type[j],".txt.rsIDmap.txt",sep=""),sep="\t")
      IDmap=IDmap[!duplicated(IDmap), ]    #remove duplicated rows
      rownames(IDmap)=IDmap[,2]
      colnames(IDmap)=c("genoposall","rsID")
      if (sum(duplicated(IDmap[,1]))>0){
        #we may have the situation that one genomic locations have multiple rsIDs (caused by 6_convertGPos2rsID.py, in Shayna's 1000G data, multiple SNPs can have the same genomic position so one sQTL can have multiple GWAS SNPs with different rsIDs)
        #we also have the situation that multiple exons can have the same SNP as their sQTL (this is not what's been tested here)
        print(paste("we have ",sum(duplicated(IDmap[,1]))," genomic locations with multiple rsIDs in ",br,sep=""))  
      }
      output=subset(output,output[,"genoposall"] %in% IDmap[,1])   #filter out sQTLs that don't have matched rsID
      output=merge(output,IDmap,all.x=TRUE,by="genoposall")  
      colnames(output)[11]="sQTL"
      #after using all.x=TRUE, if the genomic location in one row of original output matrix have, for example, 2 rsIDs, then that row will appear twice in the new output matrix and each row will be matched with one of the two rsIDs (all the rest will be the same)
      
      #sQTL rsID -> GWAS rsID
      output=merge(output,LDresult,all.x=TRUE,by="sQTL")
      
      #for each row of output matrix, check if they are related to the diseases we are interested in (we use the information in the MAPPED_TRAIT column to check the disease association)
      disease_related=rep(0,dim(output)[1])   #0 means not related and 1 means related
      for (row in 1:dim(output)[1]){
        DO=toupper(output[row,"disease.ontology"])
        if (grepl(paste(disease_key_word, collapse="|"), DO)){     #if the disease ontology contains any of the 10 disease names
          disease_related[row]=1
        }  
      }
      diseaseGWASperregion[i]=length(unique(output[which(disease_related==1),"highLD_GWAS"]))
      diseaseexonperregion[i]=length(unique(output[which(disease_related==1),"exon_full_ID"]))
      output=cbind(output,disease_related)
      
      rownames(output)=NULL
      #add exon length
      exon_length=rep(NA,dim(output)[1])
      for(k in 1:dim(output)[1]){
        exonfullinfo=strsplit(as.character(output[k,"exon_full_ID"]),split="\\|")[[1]]
        exon_length[k]=abs(as.numeric(exonfullinfo[6])-as.numeric(exonfullinfo[7]))
      }
      output=cbind(output,exon_length)
      
      #how many unique exons for each disease in each brain region?
      SNPtrait=toupper(output[,"disease.ontology"])
      for (disease in disease_key_word){
        uniqueexon=as.character(unique(output[which(grepl(disease,SNPtrait)),"exon_full_ID"]))
        diseaseexon[disease,br]=length(uniqueexon)
        numexonperdisease[[disease]]=c(numexonperdisease[[disease]],uniqueexon)   #add those unique GWAS traits in the current brain region for each disease
      }
      
      #calculate the minimum distance between SNP and splice site
      get_dis=function(outputexoninfo){
        disa=abs(as.numeric(outputexoninfo["relative_pos"])-as.numeric(strsplit(as.character(outputexoninfo["exon_full_ID"]),split="\\|")[[1]][6]))
        disb=abs(as.numeric(outputexoninfo["relative_pos"])-as.numeric(strsplit(as.character(outputexoninfo["exon_full_ID"]),split="\\|")[[1]][7]))
        return(min(disa,disb))
      }
      distance_2_SS=apply(output,1,get_dis)
      
      output=cbind(output,distance_2_SS)
      
      all_output=rbind(output,all_output)
      
      setwd(outputpath)
      write.table(output,paste("all_info_exon_sQTL_GWAS_disease_in_",br,"_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t",row.names=F)
    }
  }
  
  #fill in updated_diseaseexon and updated_diseaseGWAS
  #1.for each disease, get the unique set of GWAS loci/exon
  for (disease in disease_key_word){
    GWASloci=numGWASperdisease[[disease]][-1]   #there can be duplications (the same disease GWAS locus can show up in multiple brain regions)
    exon=numexonperdisease[[disease]][-1]       #there can be duplications (the same disease exon can show up in multiple brain regions)
    
    if (length(GWASloci)>0){
      #get the exon-SNP-GWAS loci trio based on the unique set of GWAS loci/exon 
      exon_SNPa=all_output[which(all_output[,"highLD_GWAS"] %in% GWASloci),]
      exon_SNPb=all_output[which(grepl(disease,toupper(all_output[,"disease.ontology"]))),]

      if (!identical(exon_SNPa,exon_SNPb)){   #they should be identical
        print("exon-SNP-GWAS trio not identical")
      }else{
        exon_SNP=exon_SNPb
        
        #for the current disease, store the brain region specific GWAS loci for each brain region
        disease_region_specific_GWAS=matrix(0,length(unique(unique(as.character(exon_SNP[,"highLD_GWAS"])))),length(brainregion))
        rownames(disease_region_specific_GWAS)=unique(as.character(exon_SNP[,"highLD_GWAS"]))     #this number can be greater than length(GWASloci) because when we get numGWASperdisease, for each brain region, we put the unique set of GWAS loci but those loci can come from different sQTLs
        colnames(disease_region_specific_GWAS)=brainregion
        #for the current disease, store the brain region specific exon for each brain region
        disease_region_specific_exon=matrix(0,length(unique(exon)),length(brainregion))
        rownames(disease_region_specific_exon)=unique(exon)
        colnames(disease_region_specific_exon)=brainregion
        
        #3. for each brain region, check how many of the exon-SNP pairs are significant in that brain region (the number of exon-SNP-GWAS trios may be more than the exon-SNP pairs so there can be duplicated exon-SNP pairs)
        for (BRnum in 1:length(brainregion)){
          BR=brainregion[BRnum]
          #for each exon-SNP pair, we get its p value and compare that with the permuted p value cutoff
          stored_exon=c()        #the exon that we have already checked in the current brain region (one exon can show up in multiple regions so there can be duplicated exons, when we count brain region specific exons, we don't want to double count)
          stored_GWAS=c()        #the GWAS loci that we have already checked in the current brain region (one GWAS loci can come from different sQTLs from multiple regions, any of those sQTLs is significant, the GWAS loci is considered significant, however, the count is the unique number, so we don't want to double count)
          for (pair in 1:dim(exon_SNP)[1]){
            current_exon=as.character(exon_SNP[pair,"exon_short_ID"])
            current_GWAS=as.character(exon_SNP[pair,"highLD_GWAS"])
            setwd(paste(rootoutput,"/",BR,"/","Glimmps_each_exon_cis_",BR,sep=""))
            pairresult=read.table(paste(as.character(exon_SNP[pair,"exon_short_ID"]),".asso",sep=""),sep="\t",header=TRUE)
            pairpvalue=pairresult[which(as.character(pairresult[,"SNPID"]) %in% as.character(exon_SNP[pair,"SNP_info"])),"pvals.lm"]
            
            #get the p value cutoff
            if (type[j]=="pvalue"){
              cutoff=10^-5
            }
            if (type[j]=="permutation"){
              setwd(paste(rootoutput,"/",BR,sep=""))
              cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
            }
            
            if (pairpvalue < cutoff){
              if (!(current_exon %in% stored_exon)){     #if we haven't count this exon
                updated_diseaseexon[disease,BR]=updated_diseaseexon[disease,BR]+1     
              }
              if (!(current_GWAS %in% stored_GWAS)){     #if we haven't count this GWAS loci
                updated_diseaseGWAS[disease,BR]=updated_diseaseGWAS[disease,BR]+1    
              }
              disease_region_specific_GWAS[as.character(exon_SNP[pair,"highLD_GWAS"]),BR]=1
              disease_region_specific_exon[as.character(exon_SNP[pair,"exon_full_ID"]),BR]=1
              
              stored_exon=c(stored_exon,current_exon)      #add current exon to the list (we only add significant ones)
              stored_GWAS=c(stored_GWAS,current_GWAS)      #add current GWAS locus to the list (we only add significant ones)
            }
          }
        }
      }
      setwd(outputpath)
      write.table(disease_region_specific_GWAS,paste(gsub("\\s", "_", disease),"_brain_region_specific_GWAS_loci","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
      write.table(disease_region_specific_exon,paste(gsub("\\s", "_", disease),"_brain_region_specific_exon","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
      write.table(exon_SNP,paste(gsub("\\s", "_", disease),"_exon_SNP_GWAS_trio","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
    }
  }
  
  #summary of the numbers across brain regions
  #diseaseGWAS
  unique_disease_GWAS=rep(0,length(disease_key_word))
  for (d in 1:length(disease_key_word)){
    disease=disease_key_word[d]
    unique_disease_GWAS[d]=length(unique(numGWASperdisease[[disease]]))-1
  }
  rowsum=apply(diseaseGWAS,1,sum)
  diseaseGWAS=cbind(unique_disease_GWAS,rowsum,diseaseGWAS)
  
  #updated_diseaseGWAS
  unique_disease_GWAS_updated=rep(0,length(disease_key_word))    #this should be identical to unique_disease_GWAS
  for (d in 1:length(disease_key_word)){
    disease=disease_key_word[d]
    unique_disease_GWAS_updated[d]=length(unique(numGWASperdisease[[disease]]))-1
  }
  rowsum=apply(updated_diseaseGWAS,1,sum)
  updated_diseaseGWAS=cbind(unique_disease_GWAS_updated,rowsum,updated_diseaseGWAS)
  
  #diseaseexon
  unique_disease_exon=rep(0,length(disease_key_word))
  for (d in 1:length(disease_key_word)){
    disease=disease_key_word[d]
    unique_disease_exon[d]=length(unique(numexonperdisease[[disease]]))-1
  }
  rowsum=apply(diseaseexon,1,sum)
  diseaseexon=cbind(unique_disease_exon,rowsum,diseaseexon)
  
  #unique_disease_exon_updated
  unique_disease_exon_updated=rep(0,length(disease_key_word))    #this should be identical to unique_disease_exon
  for (d in 1:length(disease_key_word)){
    disease=disease_key_word[d]
    unique_disease_exon_updated[d]=length(unique(numexonperdisease[[disease]]))-1
  }
  rowsum=apply(updated_diseaseexon,1,sum)
  updated_diseaseexon=cbind(unique_disease_exon_updated,rowsum,updated_diseaseexon)
  
  setwd(outputpath)
  write.table(diseaseGWAS,paste("GWAS_loci_in_disease_each_BR","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
  write.table(updated_diseaseGWAS,paste("updated_GWAS_loci_in_disease_each_BR","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
  write.table(diseaseexon,paste("exon_in_disease_each_BR","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
  write.table(updated_diseaseexon,paste("upadted_exon_in_disease_each_BR","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t")
  write.table(rbind(rbind(rbind(rbind(rbind(rbind(brainregion,totalsqtlnum),totalsqtlexonnum),sqtlnum),GWASnum),diseaseGWASperregion),diseaseexonperregion),paste("num_of_sQTL_and_GWASexon_and_diseaseGWASexon_per_brainregion","_",PSItype,"_",counttype,"_",splicetype,"_",type[j],".txt",sep=""),sep="\t",row.names=F,col.names=F)
}





#################################################################
#2. further summary (summarize the results above into one table)#
#################################################################
rootsqtl=strsplit(rootoutput,split=PSItype)[[1]][1]
type=c("pvalue","permutation")

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

library("biomaRt")
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'grch37.ensembl.org')
t2g=try(suppressMessages(biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name","description"), mart = mart)),silent=TRUE)   
if (!(inherits(t2g,"try-error"))){ 
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name, genetitle=description)
}else{
  print("cannot connect to BiomaRt")
  setwd("/u/home/y/yidazhan")
  t2g=as.matrix(read.table("t2g.txt",sep="\t",header=T))
}


for (t in 1:2){
  output=matrix(NA,1,28)
  colnames(output)=c("gene.symbol","exon.ID","exon.coordinate", "UCSC_Browser_link","sQTL.rsID", "GWAS.rsID", "disease.ontology", "LD.r.square",
                     paste("p.value (",formalbrainregionlist,")",sep=""), 
                     "sQTL.significant.region","significant.region","exon.length","gene.description",
                     "exon_short_ID","SNP_info","dis_2_SS")
  triolist=c()        #a list to store the exon-SNP-GWAS trait trios that we have already have the in the output matrix
  
  for (i in 1:length(brainregionlist)){
    br=brainregionlist[i]
    setwd(outputpath)
    summary=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",br,"_",PSItype,"_",counttype,"_",splicetype,"_",type[t],".txt",sep=""),sep="\t",header=T)
    
    for (j in 1:dim(summary)[1]){       #for each exon-SNP-GWAS trio
      trio=paste(summary[j,"exon_full_ID"],summary[j,"sQTL"],summary[j,"highLD_GWAS"],sep=",")
      
      if (trio %in% triolist){          #if this exon-SNP-GWAS trio is already in the table, this means that we only need to put the brain region p value and significance information into that row
        row=which(triolist %in% trio)+1   #get the row of this exon-SNP-GWAS trio (we need to add one since the first row is all NA)
        #add the current brain region into significant list (summary only contains significant result so if it is in the summary, it is significant)
        if (is.na(output[row,"sQTL.significant.region"])){       #if there is no significant brain region so far
          output[row,"sQTL.significant.region"]=formalbrainregionlist[which(brainregionlist %in% br)]
        }else{      #if we have already had significant brain regions so far
          output[row,"sQTL.significant.region"]=paste(output[row,"sQTL.significant.region"],formalbrainregionlist[which(brainregionlist %in% br)],sep=", ")
        }
        
      }else{            #if this exon-SNP-GWAS trio is not in the table, we need to add the whole thing into the table
        #get other information
        temp=strsplit(as.character(summary[j,"exon_full_ID"]),split="\\|")[[1]]
        newrow=c(temp[3],    #gene symbol
                 temp[2],    #exon Ensembl ID
                 paste(temp[4],":",temp[6],"-",temp[7],sep=""),    #exon coordinate
                 paste("=HYPERLINK(\"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",paste(temp[4],":",temp[6],"-",temp[7],sep=""),"\"",",","\"",paste(temp[4],":",temp[6],"-",temp[7],sep=""),"\")",sep=""),
                 as.character(summary[j,"sQTL"]),    #sQTL rsID  
                 as.character(summary[j,"highLD_GWAS"]),    #GWAS.rsID
                 paste(unique(strsplit(as.character(summary[j,"disease.ontology"]),split=";")[[1]]),collapse="; "),    #disease ontology (remove duplications)
                 as.character(summary[j,"r.square"]),    #LD r.square
                 rep(NA,length(brainregionlist)), #we fill in the p values later
                 #currentpvalue,    #p value
                 formalbrainregionlist[which(brainregionlist %in% br)],    #sQTL significant region (shorter name)
                 "",     #significant region
                 as.character(summary[j,"exon_length"]),    #exon length
                 strsplit(t2g[which(t2g[,"ext_gene"]==temp[3]),"genetitle"][1],split="\\[")[[1]][1],  #gene description  
                 as.character(summary[j,"exon_short_ID"]),
                 as.character(summary[j,"SNP_info"]),
                 as.numeric(summary[j,"distance_2_SS"]))   
        output=rbind(output,newrow)
        triolist=c(triolist,trio)
      }
    }
  }
  
  rownames(output)=NULL
  output=output[-1,]  #remove the first row (all NA)
  
  #fill in the p values                                       
  for (i in 1:length(brainregionlist)){   #for each brain region
    br=brainregionlist[i]
    print(br)
    for (j in 1:dim(output)[1]){   #for each unique trio
      #print(j)
      exon=output[j,"exon_short_ID"]
      snp=output[j,"SNP_info"]
      
      #get p value
      setwd(paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",br,"/","Glimmps_each_exon_cis_",br,sep=""))
      sQTLtable=read.table(paste(exon,".asso",sep=""),sep="\t",header=T)
      pvalue=sQTLtable[which(sQTLtable[,"SNPID"]==snp),"pvals.lm"]
      #get p value cutoff
      if (type[t]=="pvalue"){
        cutoff=10^-5
      }
      if (type[t]=="permutation"){
        setwd(paste(rootoutput,"/",br,sep=""))
        cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
      }
      output[j,paste("p.value (",formalbrainregionlist[which(brainregionlist %in% br)],")",sep="")]=pvalue
      if (pvalue<=cutoff){
        if (output[j,"significant.region"]==""){
          output[j,"significant.region"]=formalbrainregionlist[which(brainregionlist %in% br)]
        }else{
          output[j,"significant.region"]=paste(output[j,"significant.region"],formalbrainregionlist[which(brainregionlist %in% br)],sep=", ")
        }
      }
    }
  }
  
  setwd(outputpath)
  write.table(output,paste("sQTL_summary_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type[t],"_unsorted.txt",sep=""),sep="\t",quote=F,row.names=F)
  sortoutput=output[order(output[,1]),]   #order the output by gene symbol
  setwd(outputpath)
  write.table(sortoutput,paste("sQTL_summary_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type[t],"_sorted.txt",sep=""),sep="\t",quote=F,row.names=F)
}







