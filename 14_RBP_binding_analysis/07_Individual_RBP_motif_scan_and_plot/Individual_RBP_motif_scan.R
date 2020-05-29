#The purpose of this code:
#for RBPs with well defined motif (like NOVA and RBFOX), we perform a straight forward motif scan
#basically, we want to know that for each SNP, does it create or disrupt a binding site?
#this can also be considered as a sanity check for the deepbind model result

splicetype="SE"      #type of alternative splicing, e.g., SE, A3SS, A5SS, MXE, IR
type="pvalue"       #pvalue or permutation

outputpath=paste("/path/to/DeepBind/deepbind_individual_RBP_motif_scan",
                 splicetype,type,sep="/")
command=paste("mkdir -p",outputpath)
system(command)

#RBPs to search
RBP2scan=list("[AT]GCATG","[CT]CA[CT]")
names(RBP2scan)=c("RBFOX","NOVA")
motiflength=list(6,4)
names(motiflength)=c("RBFOX","NOVA")

#get all sQTL exons
rootinput="/path/to/summary/logit/JC"
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
sQTLexon=c()
inputpath=paste(rootinput,splicetype,sep="/")
setwd(inputpath)
for (i in 1:length(brainregionlist)){
  temp=read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregionlist[i],"_logit_JC_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)
  sQTLexon=c(sQTLexon,as.character(temp[,"exon_full_ID"]))
}
sQTLexon=unique(sQTLexon)
shortIDlist=rep(NA,length(sQTLexon))
for (e in 1:length(sQTLexon)){
  shortIDlist[e]=paste("SE",strsplit(sQTLexon[e],split="\\|")[[1]][1],sep="_")
}
exonIDconversion=cbind(shortIDlist,sQTLexon)
rownames(exonIDconversion)=shortIDlist

#read in the joblist
joblistinput=paste("/path/to/DeepBind/deepbind_input/",splicetype,"/",type,sep="")
setwd(joblistinput)
uniquejoblist=read.table(paste(splicetype,"_",type,"_uniquejoblist.txt",sep=""),sep="\t")
#get all the exon-SNP pairs
subuniquejoblist=unique(uniquejoblist[,c("Exon","SNP","Brain.region")])
rownames(subuniquejoblist)=NULL
subuniquejoblist=unique(subuniquejoblist)
rownames(subuniquejoblist)=paste(subuniquejoblist[,"Exon"],subuniquejoblist[,"SNP"],sep="~")

#generate a matrix to store the motif scan result of each pair
motifscan=matrix(NA,dim(subuniquejoblist)[1],4*length(names(RBP2scan)))
colnames(motifscan)=sort(as.vector(outer(names(RBP2scan), c("seq.wild","seq.mut","motif.num.wt","motif.num.mut"), paste, sep=".")),decreasing=T)
rownames(motifscan)=rownames(subuniquejoblist)

#get genome sequence
functionfile="/path/to/rbp_map_functions.py"
fasta="/path/to/hg19/FASTA/file"
library(reticulate)
source_python(functionfile)
genome=read_genome(fasta)
library(Biostrings)

#for each SNP, get the motif scan result
for (s in 1:dim(subuniquejoblist)[1]){
  #print(s)
  exon=as.character(subuniquejoblist[s,"Exon"])
  snp=as.character(subuniquejoblist[s,"SNP"])
  
  allele0=strsplit(snp,split="_")[[1]][3]     #reference allele, the value is 0
  allele1=strsplit(snp,split="_")[[1]][4]     #alternative allele, the value is 1
  
  if (nchar(allele0)==1 & nchar(allele1)==1){      #if it is single base mutation
    #get the position of the SNP
    chr=paste("chr",strsplit(snp,split="_")[[1]][1],sep="")
    pos=as.numeric(strsplit(snp,split="_")[[1]][2])
    strand=strsplit(sQTLexon[which(shortIDlist %in% exon)],split="\\|")[[1]][5]
    #get the sequence of the SNP
    snpseq=fetch_seq(genome,chr,pos-1,pos,strand)
    
    if (allele0==snpseq){          #the SNP is on the same strand with the exon sequence
      ref=allele0
      alt=allele1
    }else{                                        #the SNP is on the opposite strand of the exon sequence
      ref=chartr("ATGC","TACG",allele0)           #reverse complement
      alt=chartr("ATGC","TACG",allele1)
    }
    
    for (rbp in names(RBP2scan)){
      motif=as.character(RBP2scan[rbp])
      motlen=as.numeric(motiflength[rbp])
      
      #get the sequence for motif scan
      start=pos-motlen+1
      end=pos+motlen-1
      wtseq=fetch_seq(genome,chr,start-1,end,strand)
      mutseq=as.character(replaceLetterAt(DNAString(wtseq),motlen,alt))
      
      #motif scan
      motifsearchwt=gregexpr(motif, wtseq)
      motifsearchmut=gregexpr(motif, mutseq)
      
      if (length(motifsearchwt[[1]])==1 && motifsearchwt[[1]]==-1){        #there is no motif in the whole region
        wtnummotif=0
      }else{
        wtnummotif=length(motifsearchwt[[1]])      #calculate the number of matched motifs
      }
      
      if (length(motifsearchmut[[1]])==1 && motifsearchmut[[1]]==-1){        #there is no motif in the whole region
        mutnummotif=0
      }else{
        mutnummotif=length(motifsearchmut[[1]])      #calculate the number of matched motifs
      }
      
      motifscan[s,paste(rbp,"seq.wild",sep=".")]=wtseq
      motifscan[s,paste(rbp,"seq.mut",sep=".")]=mutseq
      motifscan[s,paste(rbp,"motif.num.wt",sep=".")]=wtnummotif
      motifscan[s,paste(rbp,"motif.num.mut",sep=".")]=mutnummotif
    }
  }
}
#s=232 for MAPT

fullexoninfo=exonIDconversion[as.character(subuniquejoblist[,"Exon"]),"sQTLexon"]
output=cbind(subuniquejoblist,fullexoninfo,motifscan)

setwd(outputpath)
write.table(output,paste(splicetype,type,"individual_RBP_motif_scan.txt",sep="_"),sep="\t")



#get events which create or disrupt RBFOX binding site
wtpos=which(output[,"RBFOX.motif.num.wt"]==1)
mutpos=which(output[,"RBFOX.motif.num.mut"]==1)
pos=unique(c(wtpos,mutpos))
RBFOXoutput=output[pos,]
setwd(outputpath)
write.table(RBFOXoutput,paste(splicetype,type,"individual_RBP_motif_scan_RBFOX.txt",sep="_"),sep="\t")

#get events which create or disrupt NOVA binding site
wtpos=which(output[,"NOVA.motif.num.wt"]==1)
mutpos=which(output[,"NOVA.motif.num.mut"]==1)
pos=unique(c(wtpos,mutpos))
NOVAoutput=output[pos,]
setwd(outputpath)
write.table(NOVAoutput,paste(splicetype,type,"individual_RBP_motif_scan_NOVA.txt",sep="_"),sep="\t")










