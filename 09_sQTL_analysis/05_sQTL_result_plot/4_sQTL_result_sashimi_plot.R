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
if (PSItype=="logit"){     #logit PSI with correction
  inputPSI=paste(splicetype,"_",counttype,"_logit_corrected_PSI.txt",sep="")
}
inputtotalRC=paste(splicetype,"_",counttype,"_totalRC_filter.txt",sep="")
inputIC=paste(splicetype,"_",counttype,"_IC_filter.txt",sep="")

setwd(inputpath)
PSI=read.table(inputPSI,sep="\t",header=T)              #PSI as the normalized inclusion count
setwd(totalcountinput)
totalRC=read.table(inputtotalRC,sep="\t",header=T)
IC=read.table(inputIC,sep="\t",header=T)

#transform logit PSI back to PSI
if (PSItype=="logit"){
  temp=inv.logit(as.numeric(as.matrix(PSI)))
  originalPSI=matrix(temp,dim(PSI)[1],dim(PSI)[2])
  rownames(originalPSI)=rownames(PSI)
  colnames(originalPSI)=colnames(PSI)
  PSI=as.data.frame(originalPSI)
}


#########################################
#1. single exon sQTL plot - Sashimi plot#
#########################################
#read in exon effective length table
elinput=totalcountinput
effectlengthsubset=try(suppressMessages(read.table(paste("subset_",counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t",header=T)),silent=TRUE) 
if (!(inherits(effectlengthsubset,"try-error"))){       #if we have already have the sutset file (the original file is too big)
  effectlength=effectlengthsubset
}else{      #if we don't, we need to generate that
  setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/brain_rMATS_post")
  effectlengthfull=fread(paste(counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t",header=T)
  effectlengthfull=as.data.frame(effectlengthfull)
  
  rownames(effectlengthfull)=effectlengthfull[,"ID"]
  exonID=sapply(rownames(PSI),function(x) return(strsplit(x,split="\\|")[[1]][1]))
  effectlengthsubset=effectlengthfull[exonID,]
  
  setwd(elinput)
  write.table(effectlengthsubset,paste("subset_",counttype,".raw.input.",splicetype,".txt",sep=""),sep="\t")
  effectlength=effectlengthsubset
}
rownames(effectlength)=rownames(PSI)

#read in exon annotation table
setwd("/path/to/rMATS/post/step/output")
fromGTF=read.table(paste("fromGTF.",splicetype,".txt",sep=""),sep="\t",header=T)

#generate output folder for sashimi plot
outputpathplot=paste(outputpath,"sashimi_plot",sep="/")
command=paste("mkdir -p ",outputpathplot,sep="")
system(command)


#get events and correponding information to plot
getexon=function(name){
  exoninfo=matrix(NA,length(name),11)
  rownames(exoninfo)=name
  colnames(exoninfo)=c("ID","GeneID","geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE")
  for (n in 1:length(name)){
    temp=strsplit(name[n],split="\\|")[[1]]
    
    exoninfo[n,"ID"]=temp[1]
    exoninfo[n,"GeneID"]=temp[2]
    exoninfo[n,"geneSymbol"]=temp[3]
    exoninfo[n,"chr"]=strsplit(temp[4],split="chr")[[1]][2]
    exoninfo[n,"strand"]=temp[5]
    exoninfo[n,"exonStart_0base"]=temp[6]
    exoninfo[n,"exonEnd"]=temp[7]
    exoninfo[n,"upstreamES"]=temp[8]
    exoninfo[n,"upstreamEE"]=temp[9]
    exoninfo[n,"downstreamES"]=temp[10]
    exoninfo[n,"downstreamEE"]=temp[11]
  }
  return(exoninfo)
}


for (br in 1:length(brainregionlist)){
  brainregion=brainregionlist[br]
  formalbrainregion=formalbrainregionlist[br]
  
  #1. get list of exons to plot#
  setwd(summaryinput)
  output=try(suppressMessages(read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)),silent=TRUE)   

  if (!(inherits(output,"try-error"))){       #if there are significant sQTLs
    #get the exons with disease related GWAS trait
    exontoplot=unique(as.character(output[which(output[,"disease_related"]==1),"exon_full_ID"]))
    
    #2. get the sample ID of samples in the current brain region
    setwd(paste(strsplit(rootsqtl,split="sQTL_run")[[1]],"/input_splicing/",PSItype,"/",counttype,"/",splicetype,"/",brainregion,sep=""))
    IDs.pheno=as.character(as.matrix(read.table("SRR_ID.txt",sep="\t")))
    
    #3. get the sample ID of samples in each genotype#
    if (length(exontoplot)>0){
      for (i in 1:length(exontoplot)){     #the sashimi plot needs to be done for each exon seperately
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
        
        allele0=strsplit(exongenoname,split="_")[[1]][3]     #reference allele, the value is 0
        allele1=strsplit(exongenoname,split="_")[[1]][4]     #alternative allele, the value is 1
        #0/0 -> 1, 0/1 and 1/0 -> 1, 1/1 -> 2
        Alleles<-c(paste(allele0,allele0,sep="/"),paste(allele0,allele1,sep="/"),paste(allele1,allele1,sep="/"))
        
        sampleID0=names(which(exongeno==0))
        sampleID1=names(which(exongeno==1))
        sampleID2=names(which(exongeno==2))
        samplesize=c(length(sampleID0),length(sampleID1),length(sampleID2))
        
        ###make sashimi plot###
        rootbam="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_all_bam_soft_link/"   #the folder containing all BAM files
        
        exoninfo=getexon(exon)
        #generate folders for the plot
        setwd(outputpathplot)
        outputfolder=paste(PSItype,"_",counttype,"_",splicetype,"_",type,"_",gsub("\\s", "", formalbrainregion),"_",paste(strsplit(exon,split="\\|")[[1]][-c(2,5)],collapse="_"),exongenoname,"_",exongenorsID,sep="")
        if (file.exists(outputfolder)==FALSE){   #if the folder doesn't exist, we create the folder. Otherwise, we don't do anything
          command=paste("mkdir ",outputfolder,sep="")
          system(command)
        }
        setwd(paste(outputpath,"sashimi_plot",outputfolder,sep="/"))
        
        ###prepare input files for sashimi plot###
        sortexongeno=sort(exongeno)
        sample_target=sample_control=c()  #we don't have target and control here but in order to use the previous code, we just refer to the first genotype as target and the rest as control
        sig_group=Alleles[which(samplesize!=0)][1]
        control_group=Alleles[which(samplesize!=0)][-1]
        
        grouping=matrix(NA,length(samplesize[samplesize!=0]),1)      #the file contains the grouping information of samples
        for (g in 1:length(control_group)){   
          group=control_group[g]
          #1. b1 and b2 parameter
          sample_target=names(sortexongeno[sortexongeno==which(Alleles==sig_group)-1])
          temp_control=names(sortexongeno[sortexongeno==(which(Alleles==group)-1)])
          sample_control=c(sample_control,temp_control)
          
          #2. grouping file
          grouping[1,]=paste(sig_group,": ",1,"-",length(sample_target),sep="")
          temp=as.numeric(tail(strsplit(grouping[g,],split="-")[[1]],n=1))
          grouping[g+1,]=paste(group,": ",temp+1,"-",temp+length(temp_control),sep="")
        }
        write.table(grouping,"grouping.gf",sep="\t",row.names = F,col.names = F,quote=F)
        
        b1=paste(paste(rootbam,sample_target,".bam",sep=""),collapse=",")
        b2=paste(paste(rootbam,sample_control,".bam",sep=""),collapse=",")
        
        #3. rMATS format input
        rmats=matrix(NA,dim(exoninfo)[1],23)    #the input file for sashimi plot
        colnames(rmats)=c("ID","GeneID","geneSymbol","chr","strand",
                          colnames(fromGTF)[6:11],
                          #"exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE",
                          "ID",	"IC_SAMPLE_1",	"SC_SAMPLE_1",	"IC_SAMPLE_2",	"SC_SAMPLE_2",	"IncFormLen",	"SkipFormLen",
                          "PValue",	"FDR",	"IncLevel1",	"IncLevel2",	"IncLevelDifference")
        rmats[,1:11]=exoninfo
        rmats[,"PValue"]=rmats[,"FDR"]=0
        rmats[,12]=rmats[,1]
        rownames(rmats)=rownames(exoninfo)  #we need to remove row name when outputting this into file
        
        for(e in 1:dim(rmats)[1]){
          EXON=rownames(exoninfo)[e]
          rmats[EXON,"IC_SAMPLE_1"]=paste(IC[EXON,sample_target],collapse=",")
          rmats[EXON,"SC_SAMPLE_1"]=paste((totalRC[EXON,sample_target]-IC[EXON,sample_target]),collapse=",")
          rmats[EXON,"IC_SAMPLE_2"]=paste(IC[EXON,sample_control],collapse=",")
          rmats[EXON,"SC_SAMPLE_2"]=paste((totalRC[EXON,sample_control]-IC[EXON,sample_control]),collapse=",")
          rmats[EXON,"IncFormLen"]=effectlength[EXON,"IncFormLen"]
          rmats[EXON,"SkipFormLen"]=effectlength[EXON,"SkipFormLen"]
          rmats[EXON,"IncLevel1"]=paste(PSI[EXON,sample_target],collapse=",")
          rmats[EXON,"IncLevel2"]=paste(PSI[EXON,sample_control],collapse=",")
          rmats[EXON,"IncLevelDifference"]=mean(as.matrix(PSI)[EXON,sample_target])-mean(as.matrix(PSI)[EXON,sample_control])
        }
        write.table(rmats,"rmats_events.txt",sep="\t",row.names=F,quote=F)
        
        #run sashimi plot command
        command=paste(#"source /u/local/Modules/default/init/modules.sh",
                      #"\n",
                      #"module load samtools",       
                      #"\n",
                      "/u/nobackup/yxing/PROJECT/yidazhan/research/software/anaconda2/bin/rmats2sashimiplot",
                      "--b1", b1,
                      "--b2", b2,
                      "-t", splicetype,
                      "-e", "rmats_events.txt",
                      "--l1", sig_group,
                      "--l2", "other",
                      "--exon_s", 1,
                      "--intron_s", 1,
                      "-o", "events_output",
                      "--group-info", "grouping.gf",
                      "--min-counts", 0,   
                      "--font-size", 8,
                      sep=" ")
        write.table(command,"command.sh",sep="\t",quote=F,row.names = F,col.names = F)
        command="which python"
        system(command)
        command=paste("/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -cwd -m a -M yidazhan@mail -l h_data=8G,h_rt=2:00:00 command.sh")
        system(command)
      }
    }
  }
}


