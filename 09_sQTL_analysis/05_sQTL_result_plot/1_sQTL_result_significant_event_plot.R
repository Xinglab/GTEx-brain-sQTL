args <- commandArgs(TRUE)
splicetype=args[1]
counttype=args[2]
PSItype=args[3]
inputpath=args[4]
summaryinput=args[5]
rootsqtl=args[6]
outputpath=args[7]
type=args[8]
totalcountinput=args[9]
GWASdbpath=args[10]
GWASdbname=args[11]
genoplinkprefix=args[12]
LDpath=args[13]

library("data.table")
library(gaston)
library(boot)

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

if (PSItype=="original"){   #original PSI value without any correction
  inputPSI=paste(splicetype,"_",counttype,"_PSI_filter.txt",sep="")
}
if (PSItype=="logit"){     #logit PSI with correction (because we used corrected PSI for all the sQTL calculation, when we plot the result, we also use this but only transform it back to 0-1 scale)
  inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
}
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")

setwd(inputpath)
PSI=read.table(inputPSI,sep="\t",header=T)              #PSI as the normalized inclusion count
setwd(totalcountinput)
totalRC=read.table(inputtotalRC,sep="\t",header=T)

#transform logit PSI back to PSI
if (PSItype=="logit"){
  temp=inv.logit(as.numeric(as.matrix(PSI)))
  originalPSI=matrix(temp,dim(PSI)[1],dim(PSI)[2])
  rownames(originalPSI)=rownames(PSI)
  colnames(originalPSI)=colnames(PSI)
  PSI=as.data.frame(originalPSI)
}


##########################################
#1. single exon sQTL plot - Glimmpse plot#
##########################################
print("part 1")
newoutputpath=paste(outputpath,"single_exon_Glimmpse_plot",type,sep="/")
command=paste("mkdir -p",newoutputpath)
system(command)

for (br in 1:length(brainregionlist)){
  brainregion=brainregionlist[br]
  formalbrainregion=formalbrainregionlist[br]
  
  #1. get list of exons to plot#
  setwd(summaryinput)
  output=try(suppressMessages(read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)),silent=TRUE)   
  
  if (!(inherits(output,"try-error"))){       #if there are significant sQTLs
    #get all sQTL exons
    exontoplot=unique(as.character(output[,"exon_full_ID"]))
    print(length(exontoplot))
    
    #2. get the sample ID of samples in the current brain region
    setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep=""))
    IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
    
    #3. get genotype#
    if (length(exontoplot)>0){
      for (i in 1:length(exontoplot)){
        exon=exontoplot[i]
        temp=strsplit(exon,split="\\|")[[1]]
        chrom=temp[4]      #the genotype information is by chromosome
        genotypename = paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep="")
        genoplinkfile= genoplinkprefix
        genofile = genoplinkfile
        coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
        genesymbol=temp[3]
        #.map file for the current chromosome
        c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
        map <- data.frame(c)
        names(map) <- c("Chr","SNPID","cM","Pos")
        #.raw file for the current chromosome
        d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
        genoplink <- data.frame(d)
        genoplink=genoplink[,-1] 
        IDs.geno <- names(genoplink)           #sample ID of samples with genotype information in the current brain region
        #IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]     #get the samples with genotype information
        IDs.common <- intersect(IDs.geno,IDs.pheno)
        nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
        #sub.geno <- match(IDs.common , IDs.geno)          #the position of samples in IDs.geno that are also in IDs.common
        #geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here (we basically removed the header and first column from d)
        geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs
        exongenoname=as.character(output[which(output[,"exon_full_ID"]==exon),"SNP_info"])[1]    #name of the SNP that correlates with the given exon
        exongenorsID=as.character(output[which(output[,"exon_full_ID"]==exon),"sQTL"])[1]        #rsID of the SNP that correlates with the given exon
        exongeno=geno[which(map[,"SNPID"]==exongenoname),]
        
        #4. get exon inclusion level and total read count#
        exonpsi=PSI[exon,]
        totalcount=totalRC[exon,]
        # only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
        psi=exonpsi[,IDs.common]
        count=totalcount[,IDs.common]
        
        #5. make the plot
        setwd(newoutputpath)
        # This is the same type of plot used in the Glimmps paper. 
        # It is a combination of a boxplot with scatterplot. The 
        # x-axis represents the genotypes and y-axis represents 
        # PSI values.
        allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
        allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
        #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
        Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
        
        n<- as.numeric(as.matrix(count))
        SNP<- exongeno
        Psi <- as.numeric(as.matrix(psi))
        N <- length(n)  #number of samples
        
        outfile = paste("Glimmpse_plot_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",brainregion,"_",exon,"_",exongenoname,"_",exongenorsID,".pdf",sep="")
        if (nchar(outfile)<255){     #if the file name is longer than 255 character, the file cannot be generated and also no need to generate since simple SNP won't have a name longer than 255 character
          pdf(outfile)
          
          title=paste(formalbrainregion,
                      paste("Gene: ",genesymbol,sep=""),
                      paste("Exon coordinate: ",coordinate,sep=""),
                      paste("sQTL SNP: ",exongenorsID,sep=""),
                      sep="\n")
          ylim.range <- c(0,1) #range(psi,na.rm=T)
          plot(jitter(SNP,factor=0.5), Psi,xlim=c(-0.25,2.75), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1, main=title)
          points(jitter(SNP,factor=0.5)  ,Psi  , pch= 19, cex= log10(n+1)/3 ,col=1)
          mtext(text=Alleles, side=1, at= c(0,1,2),cex=1.5,line= 1)
          points( rep(2.3,6) , seq(1,6)*0.05+0.5  , pch= 19, cex= log10( c(1,5,10,20,50,100)+1)/3  )
          text(rep(2.55,7) , seq(1,7)*0.05+0.5, c(1,5,10,20,50,100,"# reads"))
          
          par("new"=T) # add boxplot on top
          boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.75), ylim=ylim.range,at=sort(unique(SNP[!is.na(SNP)])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0),outline=FALSE)
          dev.off()
        }else{
          print(outfile)
        }
      }
    }
  }
}

#####################################################################
#2. single exon sQTL plot - Glimmpse plot (across all brain regions)#
#####################################################################
print("part 2")
newoutputpath=paste(outputpath,"Glimmpse_plot_across_all_regions",type,sep="/")
command=paste("mkdir -p",newoutputpath)
system(command)

exonIDconvert=function(shortid,fullidlist){          #given short ID, return full ID
  num_fulllist=as.character(lapply(strsplit(fullidlist,split="\\|"), `[[`, 1))
  num_shortid=strsplit(shortid,split="_")[[1]][2]
  return(fullidlist[which(num_fulllist %in% num_shortid)])
}

setwd(summaryinput)
allinfo=read.table(paste("sQTL_summary_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type,"_sorted.txt",sep=""),sep="\t",header=T,quote=NULL,check.name=F)

sQTLpair=paste(allinfo[,"exon_short_ID"],allinfo[,"SNP_info"],sep=",")       #we generate exon-SNP pair label for each row (there will be duplications)
uniquesQTLpair=unique(sQTLpair)
#get all the unique exon-SNP pairs and the corresponding information that we need
sQTLtable=matrix(NA,length(uniquesQTLpair),18)
colnames(sQTLtable)=c("exon_short_ID","exon_full_ID","sQTL.rsID","SNP_info",paste("p.value (",formalbrainregionlist,")",sep=""),"disease.ontology")


for (i in 1:dim(sQTLtable)[1]){
  sQTL=uniquesQTLpair[i]
  sQTLtable[i,"exon_short_ID"]=strsplit(sQTL,split=",")[[1]][1]
  sQTLtable[i,"exon_full_ID"]=exonIDconvert(sQTLtable[i,"exon_short_ID"],rownames(PSI))
  sQTLtable[i,"SNP_info"]=strsplit(sQTL,split=",")[[1]][2]
  sQTLtable[i,"sQTL.rsID"]=as.character(allinfo[which(allinfo[,"SNP_info"]==sQTLtable[i,"SNP_info"]),"sQTL.rsID"])[1]
  sQTLtable[i,paste("p.value (",formalbrainregionlist,")",sep="")]=as.numeric(allinfo[which(sQTLpair==sQTL)[1],paste("p.value (",formalbrainregionlist,")",sep="")])
  sQTLtable[i,"disease.ontology"]=paste(unique(unlist(strsplit(as.character(allinfo[which(sQTLpair==sQTL),"disease.ontology"]),split="; "))),collapse=";")  
}
print(dim(sQTLtable)[1])

for (i in 1:dim(sQTLtable)[1]){   #for each sQTL
  print(i)
  
  exon=sQTLtable[i,"exon_full_ID"]
  snprsid=as.character(sQTLtable[i,"sQTL.rsID"])
  snpinfo=as.character(sQTLtable[i,"SNP_info"])
  GWAStrait=as.character(sQTLtable[i,"disease.ontology"])
  
  setwd(newoutputpath)
  outfile=paste("Glimmpse_plot_all_sQTL_all_brainregion_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",exon,"_",snprsid,"_",snpinfo,"_",i,".pdf",sep="")
  if (nchar(outfile)<255){
    pdf(outfile)
    
    par(mfrow = c(3, 5)) 
    par(mar = c(3, 3.5, 1, 0.5))   #margin of each individual plot
    par(mgp = c(1.5, 0.5, 0))
    par(oma = c(0, 0, 7, 1))   #outside margin (keep some room for the main title)
    
    for (br in 1:length(brainregionlist)){
      brainregion=brainregionlist[br]
      formalbrainregion=formalbrainregionlist[br]
      pvalue=as.character(sQTLtable[i,paste("p.value (",formalbrainregion,")",sep="")])
      
      #get the sample ID of samples in the current brain region
      setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep=""))
      IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
      
      #get exon information as title
      temp=strsplit(exon,split="\\|")[[1]]
      chrom=temp[4]      #the genotype information is by chromosome
      genotypename = paste(rootsqtl,"/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep="")
      genoplinkfile= genoplinkprefix
      genofile = genoplinkfile
      coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
      genesymbol=temp[3]
      genetitle=allinfo[which(allinfo[,"gene.symbol"]==genesymbol),"gene.description"][1]
      
      #.map file for the current chromosome
      c <-fread(paste(genotypename,"/",chrom,".",genofile, ".map", sep=""), header=F, sep="\t")    # Note: check sep (My code)
      map <- data.frame(c)
      names(map) <- c("Chr","SNPID","cM","Pos")
      #.raw file for the current chromosome
      d <- fread(paste(genotypename,"/",chrom, ".", genofile, ".map_tpose.raw", sep=""), header=T,sep="\t",  na.strings=c("NA"))     # (My code)
      genoplink <- data.frame(d)
      genoplink=genoplink[,-1] 
      IDs.geno <- names(genoplink)           #sample ID of samples with genotype information in the current brain region
      #IDs.common <- IDs.geno[is.element(IDs.geno, IDs.pheno)]     #get the samples with genotype information
      IDs.common <- intersect(IDs.geno,IDs.pheno)
      nsnps <- dim(genoplink)[1] -5                              #remove the first a few lines
      #sub.geno <- match(IDs.common , IDs.geno)          #the position of samples in IDs.geno that are also in IDs.common
      #geno <- as.matrix(genoplink[seq(1, nsnps)+5, sub.geno])    #  I changed here (we basically removed the header and first column from d)
      geno <- as.matrix(genoplink[seq(1, nsnps)+5, IDs.common])   #this is the genotype table for all SNPs
      
      #get genotype information
      exongenoname=snpinfo    #name of the SNP that correlates with the given exon
      exongenorsID=snprsid        #rsID of the SNP that correlates with the given exon
      exongeno=geno[which(map[,"SNPID"]==exongenoname),]
      
      #get exon inclusion level and total read count#
      exonpsi=PSI[exon,]
      totalcount=totalRC[exon,]
      # only take those individuals with genotypes from the count table, and sort the phenotype to be the same order as in the genotype matrix
      psi=exonpsi[,IDs.common]
      count=totalcount[,IDs.common]
      
      #5. make the plot
      # This is the same type of plot used in the Glimmps paper. 
      # It is a combination of a boxplot with scatterplot. The 
      # x-axis represents the genotypes and y-axis represents 
      # PSI values.
      allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
      allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
      #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
      Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
      n<- as.numeric(as.matrix(count))
      SNP<- exongeno
      Psi <- as.numeric(as.matrix(psi))
      N <- length(n)  #number of samples
      
      title=paste(formalbrainregion,
                  paste("P value: ",pvalue,sep=""),
                  sep="\n")
      
      #get the p value cutoff
      if (type=="pvalue"){
        cutoff=10^-5
      }
      if (type=="permutation"){
        setwd(paste(rootsqtl,PSItype,counttype,splicetype,brainregion,sep="/"))
        cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
      }
      if (as.numeric(pvalue)<=cutoff){
        color.title="red"
      }else{
        color.title="black"
      }
      
      ylim.range <- c(0,1) #range(psi,na.rm=T)
      plot(jitter(SNP,factor=0.5), Psi,xlim=c(-0.25,2.75), ylab="",xlab="",xaxt="n",type="n" ,ylim=ylim.range, cex.main=1)
      points(jitter(SNP,factor=0.5)  ,Psi  , pch= 19, cex= log10(n+1)/1 ,col=1)
      mtext(text=title, side=3, line=0.2, col=color.title, cex=0.7)
      mtext(text=Alleles, side=1, at= c(0,1,2),cex=0.7,line= 0.3)
      #points( rep(2.3,6) , seq(1,6)*0.05+0.5  , pch= 19, cex= log10( c(1,5,10,20,50,100)+1)/1  )
      #text(rep(2.55,7) , seq(1,7)*0.05+0.5, c(1,5,10,20,50,100,"# reads"))
      par("new"=T) # add boxplot on top
      boxplot(Psi~SNP, ylab="", xlab="", xaxt="n", yaxt="n",boxwex=0.35, xlim=c(-0.25,2.75), ylim=ylim.range,at=sort(unique(SNP[!is.na(Psi)][!is.na(SNP[!is.na(Psi)])])) ,border=gray(0.45),col=rgb(1,1,1,alpha=0.6),outline=FALSE)  
    }
    main.title1=paste(genesymbol,coordinate,exongenorsID,sep="      ")
    main.title2=genetitle
    main.title3=GWAStrait
    if (is.na(main.title2)){
      main.title=paste(main.title1,main.title3, sep="\n")
    }else{
      main.title=paste(main.title1,main.title2,main.title3, sep="\n")
    }
    mtext(main.title, outer = TRUE,side = 3, cex = 0.7, line = 2)
    dev.off()
  }else{
    print(outfile)
  }
}



###########################################
#3. single exon sQTL plot - Haploview plot#
###########################################
run=FALSE       #don't run this part
if (run){
  setwd(GWASdbpath)
  GWAS_mine=read.table(GWASdbname,sep="\t",header=T,fill = TRUE,quote=NULL)
  
  disease_key_word=toupper(c("Alzheimer","Amyotrophic lateral sclerosis","Parkinson","frontotemporal dementia","Huntington",
                             "epilepsy","autism","schizophrenia","bipolar","depression",
                             "attention deficit hyperactivity disorder","glio","multiple sclerosis","narcolepsy","stroke"))       #change all the names to upper case for comparison
  disease_short_kw=c("AD","ALS","PD","FTD","HD","Epilepsy","Autism","Schizophrenia","Bipolar","Depression",
                     "ADHD","Glioma/Glioblastoma","MS","Narcolepsy","Stroke")
  
  rs2disease=function(rsid){         #get the disease name of the given SNP (this only applies for SNPs that are related to at least one of the 10 diseases)
    disease_name=rep(NA,length(rsid))
    for (i in 1:length(rsid)){   #for each rsID
      ID=rsid[i]
      #get the disease trait information of the current SNP
      annotation=as.character(output[which(output[,"highLD_GWAS"]==ID),"disease.ontology"][1])
      #check the diseases that are in the annotation
      diseaselist=c()
      for (j in 1:length(disease_key_word)){
        if (grepl(disease_key_word[j],toupper(annotation))){
          diseaselist=c(diseaselist,disease_short_kw[j]) 
        }
      }
      disease_name[i]=paste("(",paste(diseaselist,collapse=","),")",sep="")
    }
    return(disease_name)
  }
  
  for (br in 1:length(brainregionlist)){
    brainregion=brainregionlist[br]
    formalbrainregion=formalbrainregionlist[br]
    
    #1. get list of exons to plot#
    setwd(summaryinput)
    output=try(suppressMessages(read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)),silent=TRUE)   
    #output=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
    
    if (!(inherits(output,"try-error"))){       #if there are significant sQTLs
      #get the exons with disease related GWAS trait
      exontoplot=unique(as.character(output[which(output[,"disease_related"]==1),"exon_full_ID"]))
      
      setwd(outputpath)
      if (length(exontoplot)>0){
        for (i in 1:length(exontoplot)){      #for each exon
          exon=exontoplot[i]
          chr=strsplit(strsplit(exon,split="\\|")[[1]][4],split="chr")[[1]][2]
          temp=strsplit(exon,split="\\|")[[1]]
          coordinate=paste(temp[4],":",temp[6],"-",temp[7],sep="")
          genesymbol=temp[3]
          
          #get the rsID of sQTL SNP and rsID of all the other SNPs that are in highLD with it
          temp=output[which(output[,"exon_full_ID"]==exon),]
          rsID=c(unique(as.character(temp[,"sQTL"])),unique(as.character(temp[,"highLD_GWAS"])))          #rsID of sQTL SNP + rsID of all the other SNPs that are in highLD with it
          if (length(rsID)==2 && rsID[1]==rsID[2]){       #if the sQTL is in high LD with itself and there is no other GWAS loci
            donothing="donothing"     #we don't do anything in this situation
          }else{
            disease_related=temp[,"disease_related"]      #disease status of those GWAS SNPs
            label=c("(sQTL)",disease_related)        #label of each SNP
            label[which(label==1)]=rs2disease(rsID[which(label==1)])
            #label[which(label==1)]="(disease)"
            label[which(label==0)]=" "
            
            #remove SNPs that are not related to any diseases
            rsID=rsID[which(label!=" ")]
            label=label[which(label!=" ")]
            
            setwd(LDpath)
            #LDpair=read.table(paste("LD_chr",chr,"_within_window.ld",sep=""),header=T)
            LDpair=fread(paste("LD_chr",chr,"_within_window.ld",sep=""),header=T)
            
            #get genomic position and r square
            pos=rep(NA,length(rsID))
            for (j in 1:length(rsID)){ #for each snp
              snp=rsID[j]
              tempLDpair=LDpair    #we generate a copy of the original data
              #let's look at SNP_A and see if snp is in the list
              setkey(tempLDpair, SNP_A)
              subsettemp=tempLDpair[.(snp)]      #look for rows in tempLDpair that have the SNP_A value equal the current snp
              if (dim(subsettemp)[1]>0 && as.character(unique(subsettemp[,"BP_A"]))!="NA"){    #if our snp is in SNP_A   
                #of note, tempLDpair[.(snp)] will return one row no matter the snp is in SNP_A or not so dim(subsettemp)[1]>0 will always be true
                #but if as.character(unique(subsettemp[,"BP_A"]))!="NA", it means that the function only returns one row and that row doesn't contain any information, i.e., snp is not in SNP_A
                pos[j]=as.character(unique(subsettemp[,"BP_A"]))
              }else{           #if our snp is not in SNP_A
                tempLDpair=LDpair    #we generate a copy of the original data
                #let's look at SNP_B and see if snp is in the list
                setkey(tempLDpair, SNP_B)
                subsettemp=tempLDpair[.(snp)]
                if (dim(subsettemp)[1]>0){        
                  #like what is said above, this could still give NA but here it is fine because we have already searched SNP_A
                  #if there is no match in SNP_B, then there is no other place to search so if it is NA, then it means that the snp is not in the list
                  pos[j]=as.character(unique(subsettemp[,"BP_B"]))
                }
              }
            }
            pos[which(pos=="NA")]=0       #if the position is still missing, we just use 0
            #rank the rsID/pos/label based on position
            rsID=rsID[order(pos)]
            label=label[order(pos)]
            pos=pos[order(pos)]
            
            #get pairwise r-square and put it into a square matrix
            haploview=matrix(NA,length(rsID),length(rsID))
            rownames(haploview)=colnames(haploview)=paste(rsID,label,sep=" ")
            for (a in 1:(length(rsID)-1)){      #we only fill in the upper triangle matrix (because the diagnal is all 0 and the lower triangle matrix is the same)
              #for (a in 1:length(rsID)){
              for (b in (a+1):length(rsID)){
                #for (b in 1:length(rsID)){
                snpa=rsID[a]
                snpb=rsID[b]
                
                if (snpa==snpb){       #this could happen, i.e., a sQTL is in high LD with itself (and this is not true for all sQTLs)
                  haploview[a,b]=1
                }else{
                  #situation 1: snpa is in SNP_A and snpb is in SNP_B
                  tempLDpair=LDpair
                  setkey(tempLDpair, SNP_A)
                  subsettemp=tempLDpair[.(snpa)]     #search snpa in SNP_A
                  if (dim(subsettemp)[1]>0 && as.character(unique(subsettemp[,"R2"]))!="NA"){          #if snpa is in SNP_A
                    if (snpb %in% as.matrix(subsettemp[,"SNP_B"])){          #if snpb is in SNP_B
                      haploview[a,b]=as.numeric(subsettemp[which(subsettemp[,"SNP_B"]==snpb),"R2"])
                    }
                  }
                  
                  #situation 2: snpa is in SNP_B and snpb is in SNP_A
                  if (is.na(haploview[a,b])){         #if the first situation is not true
                    tempLDpair=LDpair
                    setkey(tempLDpair, SNP_B)
                    subsettemp=tempLDpair[.(snpa)]     #search snpa in SNP_B
                    if (dim(subsettemp)[1]>0 && as.character(unique(subsettemp[,"R2"]))!="NA"){          #if snpa is in SNP_B
                      if (snpb %in% as.matrix(subsettemp[,"SNP_A"])){          #if snpb is in SNP_A
                        haploview[a,b]=as.numeric(subsettemp[which(subsettemp[,"SNP_A"]==snpb),"R2"])
                      }
                    }
                  }
                  
                }
              }
            }
            
            #fill in other part of the matrix
            diag(haploview)=1
            for (a in 2:length(rsID)){     
              for (b in 1:(a-1)){
                haploview[a,b]=haploview[b,a]
              }
            }
            
            haploview[is.na(haploview)]=0   #if LD r-square is missing, we use 0 (missing means r-square < 0.8)
            
            #make the plot (we use LD.plot in gaston package instead of LDheatmap in LDheatmap package)
            #LD.plot: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/LD.plot
            
            setwd(outputpath)
            outfile = paste("Haploview_plot_",PSItype,"_",counttype,"_",splicetype,"_",type,"_",brainregion,"_",exon,"_",rsID[which(label=="(sQTL)")],".pdf",sep="")
            pdf(outfile,width=8,height=10)
            
            par(mai=c(0,0.82,0,0.42))
            title=paste(paste(formalbrainregion,
                              paste("Gene: ",genesymbol,sep=""),
                              paste("sQTL SNP: ",rsID[which(label=="(sQTL)")],sep=""),
                              sep="    "),
                        paste("Exon coordinate: ",coordinate,sep=""),
                        sep="\n")
            LD.plot(haploview,snp.positions=as.numeric(pos),
                    graphical.par = list(cex = 0.5, mar = c(0, 0, 0, 0), col="red"),   
                    #text size and color (the text on each cell and the text for SNP label/title of plot must share the same color so we just use red here
                    #because white or black can only work for one place but not both (the background for each cell is black but the background for the plot is white).
                    #since we need to manually edit the plot anyway, we just use red for now and we can change the color later manually)
                    #cex.ld = 4,
                    cex.snp = 1.5,
                    polygon.par = list(border = "white"),    #add white border for each cell
                    draw.chr=TRUE,
                    color.scheme = function(ld) rgb(1-abs(ld),1-abs(ld),1-abs(ld),maxColorValue=1))   
            #the function here specify a range of value from rgb(1,1,1) when ld=0, which is white to rgb(0,0,0) when ld=1, which is black
            #ld is the number in the haploview matrix
            mtext(text=title, side=1,cex=1,line=-1.5)
            dev.off()
          }
        }
      }
    }
  }
}

print("finished!")
