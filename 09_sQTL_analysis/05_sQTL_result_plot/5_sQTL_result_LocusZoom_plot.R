#to do: add BED file result into the plot

args <- commandArgs(TRUE)
splicetype=args[1]
counttype=args[2]
PSItype=args[3]
inputpath=args[4]
summaryinput=args[5]
rootsqtl=args[6]
outputpath=args[7]
type=args[8]
IDconvertscript=args[9]
locuszoom=args[10]
extendwide=30000      #we extend the plot on each side to include the information that falls into the two boundaries
extendnarrow=30

library("data.table")
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

checkpos=function(SNPpos,exon,splicetype){
  temp=strsplit(exon,split="\\|")[[1]]
  strand=temp[5]
  if (splicetype=="SE"){
    exonstart=as.numeric(temp[6])
    exonend=as.numeric(temp[7])
    if (strand=="+"){
      SS5start=exonend-3
      SS5end=exonend+6
      SS3start=exonstart-20
      SS3end=exonstart+3
    }
    if (strand=="-"){
      SS5start=exonstart-6
      SS5end=exonstart+3
      SS3start=exonend-3
      SS3end=exonend+20
    }
  }
  if (splicetype=="SE"){
    if (SNPpos>=exonstart && SNPpos<=exonend){       #if the SNP is on exon body, the distance is 0
      dis=0
    }else{
      dis=min(abs(SNPpos-SS5start),abs(SNPpos-SS5end),abs(SNPpos-SS3start),abs(SNPpos-SS3end))  
    }
    if (SNPpos>=SS5start && SNPpos<=SS5end){
      label="SS"
    }else if (SNPpos>=SS3start && SNPpos<=SS3end){
      label="SS"
    }else {
      label="not_SS"
    }
  }
  
  if (splicetype=="A5SS"){
    exonstart=as.numeric(temp[6])     #we choose the longer exon
    exonend=as.numeric(temp[7])
    if (strand=="+"){
      ass1=as.numeric(temp[7])
      ass1_start=ass1-3
      ass1_end=ass1+6
      ass2=as.numeric(temp[9])
      ass2_start=ass2-3
      ass2_end=ass2+6
    }
    if (strand=="-"){
      ass1=as.numeric(temp[6])
      ass1_start=ass1-6
      ass1_end=ass1+3
      ass2=as.numeric(temp[8])
      ass2_start=ass2-6
      ass2_end=ass2+3
    }
  }
  if (splicetype=="A3SS"){
    exonstart=as.numeric(temp[6])     #we choose the longer exon
    exonend=as.numeric(temp[7])
    if (strand=="+"){
      ass1=as.numeric(temp[6])
      ass1_start=ass1-20
      ass1_end=ass1+3
      ass2=as.numeric(temp[8])
      ass2_start=ass2-20
      ass2_end=ass2+3
    }
    if (strand=="-"){
      ass1=as.numeric(temp[7])
      ass1_start=ass1-3
      ass1_end=ass1+20
      ass2=as.numeric(temp[9])
      ass2_start=ass2-3
      ass2_end=ass2+20
    }
  }   
  if (splicetype=="A3SS" || splicetype=="A5SS"){
    if (SNPpos>=exonstart && SNPpos<=exonend){       #if the SNP is on exon body, the distance is 0
      dis=0
    }else{
      dis=min(abs(SNPpos-ass1_start),abs(SNPpos-ass1_end),abs(SNPpos-ass2_start),abs(SNPpos-ass2_end)) 
    }
    if (SNPpos>=ass1_start && SNPpos<=ass1_end){
      label="SS"
    }else if (SNPpos>=ass2_start && SNPpos<=ass2_end){
      label="SS"
    }else {
      label="not_SS"
    }
  }
  
  return(list(distance=dis,label=label))
}


####################################
#1. single sQTL exon locuszoom plot#
####################################
newoutputpath=paste(outputpath,"locus_zoom_plot",type,sep="/")
command=paste("mkdir -p",newoutputpath)
system(command)

for (br in 1:length(brainregionlist)){
  brainregion=brainregionlist[br]
  formalbrainregion=formalbrainregionlist[br]
  
  #1. get list of exons to plot#
  setwd(summaryinput)
  output=try(suppressMessages(read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)),silent=TRUE)   
  output_rowname=paste(as.character(output[,"exon_full_ID"]),as.character(output[,"sQTL"]),sep="~")
  
  if (!(inherits(output,"try-error"))){       #if there are significant sQTLs
    exonSNPtoplot=unique(paste(as.character(output[,"exon_full_ID"]),as.character(output[,"sQTL"]),sep="~"))     #unique exon pair
    exonlist=unlist(sapply(exonSNPtoplot,strsplit,split="~"))[seq(1,length(exonSNPtoplot)*2-1,2)]                #corresponding exon ID
    print(length(exonSNPtoplot))
    
    if (length(exonSNPtoplot)>0){
      for (i in 1:length(exonSNPtoplot)){
        exon=exonlist[i]
        temp=strsplit(exon,split="\\|")[[1]]
        chr=temp[4]
        genesymbol=temp[3]
        exonshortID=paste("SE",temp[1],sep="_")
        #get all the SNPs within 200kb of the current exon & the corresponding p value of those SNPs
        path=paste(rootsqtl,PSItype,counttype,splicetype,brainregion,paste("Glimmps_each_exon_cis_",brainregion,sep=""),sep="/")
        setwd(path)
        cortable=read.table(paste(exonshortID,".asso",sep=""),sep="\t",header=T)
        
        #for each SNP, check whether it is on the splice sites of the current exon + calculate the distance of that to the nearest splice site
        SS=dis=rep(NA,dim(cortable)[1])
        for (s in 1:dim(cortable)[1]){
          result=checkpos(cortable[s,"Pos"],exon,splicetype)
          SS[s]=result$label
          dis[s]=result$distance
        }
        cortable=cbind(cortable,SS,dis)
        
        if (type=="pvalue"){
          cutoff=10^-5
        }
        if (type=="permutation"){
          setwd(paste(rootsqtl,PSItype,counttype,splicetype,brainregion,sep="/"))
          cutoff=as.numeric(as.matrix(read.table("permutation_p_value_cutoff.txt",sep="\t")))
        }
        
        for (size in c("wide","narrow")){
          if (size=="wide"){          #we plot everything within 200kb
            width=Inf
            extend=extendwide
          }
          if (size=="narrow"){        #we only plot SNPs within 300b
            width=300
            extend=extendnarrow
          }
          
          #1. get all the SNP ID (get rsID from genopos)
          subcortable=subset(cortable,cortable[,"dis"]<=width)
          if (dim(subcortable)[1]>0){   #if we have SNPs within the range
            snpconvert=matrix(NA,dim(subcortable)[1],5)
            rowname=rep(NA,dim(subcortable)[1])
            pos=rep(NA,dim(subcortable)[1])
            for (s in 1:dim(subcortable)[1]){
              snpconvert[s,1]=exonshortID
              snptemp=strsplit(as.matrix(subcortable[s,"SNPID"]),split="_")[[1]]
              snpconvert[s,2]=subcortable[s,"Chr"]
              snpconvert[s,3]=as.character(as.matrix(subcortable[s,"SNPID"]))
              snpconvert[s,4]=as.numeric(snptemp[2])
              snpconvert[s,5]=subcortable[s,"pvals.lm"]
              rowname[s]=paste(snptemp[1],snptemp[2],sep="_")
              pos[s]=as.numeric(snptemp[2])
            }
            sigsubcortable=subset(subcortable,subcortable[,"pvals.lm"]<=cutoff)
            sigSSsubcortable=subset(sigsubcortable,sigsubcortable[,"SS"]=="SS")
            if (dim(sigSSsubcortable)[1]>0){   #if we have significant SNPs on splice site
              opath=paste(newoutputpath,paste(paste("locuszoom_SpliceSite_",size,"_",type,"_",brainregion,sep=""),gsub("\\|",",",exonSNPtoplot[i]),sep="~"),sep="/")
            }else{
              opath=paste(newoutputpath,paste(paste("locuszoom_",size,"_",type,"_",brainregion,sep=""),gsub("\\|",",",exonSNPtoplot[i]),sep="~"),sep="/")
            }
            
            command=paste("mkdir -p",opath)
            system(command)
            setwd(opath)
            
            setwd(opath)
            write.table(snpconvert,"snpconvert.txt",sep="\t",quote=F,row.names=F,col.names=F)
            
            #2. convert them into rsID
            command=paste("/path/to/python",IDconvertscript,paste(opath,"snpconvert.txt",sep="/"),strsplit(chr,split="chr")[[1]][2])
            system(command)
            
            #3. generate input for locus zoom plot
            subcortable=as.matrix(subcortable)
            rownames(subcortable)=rowname
            rsID=rep(NA,dim(subcortable)[1])
            IDconversion=try(suppressMessages(as.matrix(read.table("snpconvert.txt.rsIDmap.txt",sep="\t"))),silent=TRUE)
            if (!(inherits(IDconversion,"try-error"))){          #if we have SNPs after the ID conversion (sometimes we will have no SNP with rsID)
              lz=matrix(,0,2)
              for (s in 1:dim(subcortable)[1]){       #for each SNP within the range, if we can find its rsID, we add it to the input file
                if (rowname[s] %in% IDconversion[,1]){
                  rsID[s]=IDconversion[which(IDconversion[,1] %in% rowname[s])[1],2]
                  newrow=c(IDconversion[which(IDconversion[,1] %in% rowname[s])[1],2],subcortable[s,"pvals.lm"])
                  lz=rbind(lz,newrow)
                }else{
                  rsID[s]="ID_not_found"
                }
              }
              subcortable=cbind(subcortable,rsID)
              colnames(lz)=c("MarkerName","P-value")
              setwd(opath)
              write.table(lz,"locuszoom_input.txt",sep="\t",quote=F,row.names=F)
              
              #4. generate SNP label file (sQTL SNP + significant SNPs on splice site)
              label=matrix(,0,3)
              colnames(label)=c("snp",	"string",	"color")
              #add sQTL SNP if it is within the range
              sQTLrsID=as.character(as.matrix(output[which(output_rowname %in% exonSNPtoplot[i]),"sQTL"]))[1]  #rsID of the sQTL SNP
              if (sQTLrsID %in% lz[,1]){       #if the sQTL SNP is within the plot range
                #if (subcortable[IDconversion[which(IDconversion[,2] %in% as.character(as.matrix(output[which(output_rowname %in% exonSNPtoplot[i]),"sQTL"])))[1],1],"SS"]=="SS"){
                if (subcortable[which(subcortable[,"rsID"]==sQTLrsID),"SS"]=="SS"){  
                  #if the sQTL is also on SS
                  newrow=c(as.character(as.matrix(output[which(output_rowname %in% exonSNPtoplot[i]),"sQTL"]))[1],"sQTL-SS","blue")
                  label=rbind(label,newrow)
                }else{
                  newrow=c(as.character(as.matrix(output[which(output_rowname %in% exonSNPtoplot[i]),"sQTL"]))[1],"sQTL","blue")
                  label=rbind(label,newrow)
                }   
              }
              #add all other significant SNPs on SS
              sigsubcortable=subset(subcortable,as.numeric(subcortable[,"pvals.lm"])<=cutoff)
              sigSSsubcortable=subset(sigsubcortable,sigsubcortable[,"SS"]=="SS")
              if (dim(sigSSsubcortable)[1]>0){      #if we have significant SNPs on SS
                for (sigSS in 1:dim(sigSSsubcortable)[1]){     #add each of them into the label matrix
                  if (dim(label)[1]>0 && (sigSSsubcortable[sigSS,"rsID"] %in% label[1,"snp"])){       #if this SNP is already in the list (also a sQTL SNP)
                    next
                  }else{
                    newrow=c(as.matrix(sigSSsubcortable[sigSS,"rsID"]),"SS","purple")
                    label=rbind(label,newrow)
                  }
                }
              }          
              setwd(opath)
              write.table(label,"locuszoom_SNP_label.txt",sep="\t",quote=F,row.names=F)
              
              
              #5. run locuszoom
              setwd(opath)
              if (dim(label)[1]>0){         #if we have anything we want to label
                command=paste(locuszoom,"--metal locuszoom_input.txt --denote-markers-file locuszoom_SNP_label.txt",
                              "--pop EUR --build hg19 --source 1000G_March2012",
                              "--chr",strsplit(chr,split="chr")[[1]][2],"--start",range(pos)[1]-extend,"--end",range(pos)[2]+extend,
                              #"--gwas-cat whole-cat_significant-only",
                              "--plotonly",
                              paste("signifLine=",-log10(cutoff),sep=""),
                              "signifLineColor=blue")
              }else{                        #if we don't have anything we want to label
                command=paste(locuszoom,"--metal locuszoom_input.txt",
                              "--pop EUR --build hg19 --source 1000G_March2012",
                              "--chr",strsplit(chr,split="chr")[[1]][2],"--start",range(pos)[1]-extend,"--end",range(pos)[2]+extend,
                              #"--gwas-cat whole-cat_significant-only",
                              "--plotonly",
                              paste("signifLine=",-log10(cutoff),sep=""),
                              "signifLineColor=blue")
              }
              system(command)
              #6. output all the information 
              setwd(opath)
              write.table(subcortable,"subcortable.txt",sep="\t",quote=F)
              write.table(output[output_rowname %in% exonSNPtoplot[i],],"suboutput.txt",sep="\t")
            }
          }
        }
      }
    }
  }
}
