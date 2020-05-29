setwd("/path/to/sQTL/scripts")
brainregion=as.character(as.matrix(read.table("brainregionlist.txt",sep="\t")))
inputpath="/output/path/of/1_Fraction_of_exon_within_300bp"
outputpath=outputpath


splicetypelist=c("SE","A3SS","A5SS")
for (splicetype in splicetypelist){
  #make the plot
  library(ggplot2)
  setwd(inputpath)
  proportionmatrix=matrix(,nrow=0,ncol=3)
  for (i in 1:length(brainregion)[1]){
    br=brainregion[i]
    proportionresult=as.matrix(read.table(paste(splicetype,"_",br,"_proportion.txt",sep=""),sep="\t",header=T))
    brain.region=rep(br,dim(proportionresult)[1])
    temp=cbind(proportionresult[,c("proportion","xaxis")],brain.region)
    proportionmatrix=rbind(proportionmatrix,temp)
  }
  colnames(proportionmatrix)=c("proportion","transformed.p.value.cutoff","brain.region")
  rownames(proportionmatrix)=NULL
  
  #change brain region labels
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-CerebellarHemisphere"),"brain.region"]="Cerebellar Hemisphere"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Cerebellum"),"brain.region"]="Cerebellum"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Cortex"),"brain.region"]="Cortex"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-FrontalCortexBA9"),"brain.region"]="Frontal Cortex"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-AnteriorcingulatecortexBA24"),"brain.region"]="Anterior cingulate cortex"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Hypothalamus"),"brain.region"]="Hypothalamus"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Caudatebasalganglia"),"brain.region"]="Caudate"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Amygdala"),"brain.region"]="Amygdala"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Spinalcordcervicalc-1"),"brain.region"]="Spinal cord"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Nucleusaccumbensbasalganglia"),"brain.region"]="Nucleus accumbens"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Substantianigra"),"brain.region"]="Substantia nigra"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Putamenbasalganglia"),"brain.region"]="Putamen"
  proportionmatrix[which(proportionmatrix[,"brain.region"]=="Brain-Hippocampus"),"brain.region"]="Hippocampus"
  
  df <- as.data.frame(proportionmatrix)
  df$transformed.p.value.cutoff=as.numeric(as.matrix(df$transformed.p.value.cutoff))
  df$proportion=as.numeric(as.matrix(df$proportion))
  start=range(df[,"transformed.p.value.cutoff"])[1]
  end=range(df[,"transformed.p.value.cutoff"])[2]
  
  setwd(outputpath)
  pdf(paste(splicetype,"fraction of exons with sQTL SNPs<300bp.pdf"),height=8,width=8)
  cbPalette <- c("#E53537",     #GTEx Artery color/Amygdala
                 "#EEA760",     #GTEx adipose tissue color/ACC
                 "#C9E5C3",     #GTEx pituitary color/Caudate
                 "#9ACB3C",     #GTEx lung color/CH
                 "#028B45",     #GTEx thyroid color/Cerebellum
                 "#FFD923",     #GTEx nerve color/cortex
                 "#F0EC68",     #GTEx brain color/frontal cortex
                 "#6D67AF",     #GTEx muscle color/hippocampus
                 "#7F388D",     #GTEx Heart color/Hypothalamus
                 "#2999D5",     #GTEx Skin color/Nucleus accumbens
                 "#4EC3C7",     #GTEx breast color/Putamen
                 "#CD8ABC",     #GTEx LCL color/spinal cord
                 "#F287B7")     #GTEx blood color/Substantia nigra
  p=ggplot(data=df, aes(x=transformed.p.value.cutoff, y=proportion, group=brain.region)) +
    geom_line(aes(color=brain.region),size = 1.5)+
    scale_color_manual(values=cbPalette) +
    labs(x=expression('-log'[10]*'(P)'), y = "Fraction of exons with sQTL SNPs<300bp") + 
    #scale_x_continuous(limits = c(floor(quantile(-log10(allpvalue),c(0.01,0.99))[1]),ceiling(quantile(-log10(allpvalue),c(0.01,0.99))[2]))) + 
    scale_x_continuous(breaks=seq(2,end,by=2),limits = c(2,end)) +       #we start from 2 instead of 1 (this will create warnings because of rows with transformed.p.value.cutoff<2. We don't need to worry about that)
    scale_y_continuous(limits=c(0,1)) + 
    theme(
      # axis
      axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
      axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
      axis.text = element_text(size = rel(1.1)),
      axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
      axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
      axis.line = element_line(colour = "black"),
      
      #background
      panel.background = element_blank(),
      #panel.grid.minor = element_line(colour = "grey"),
      
      #legend
      legend.title = element_text(face = "bold"),
      legend.key = element_rect(fill = "white", colour = "black"),
      legend.text = element_text(size = 10, face = "bold"),
      
      # strip
      strip.text=element_text(size = rel(1.3)),
      aspect.ratio=1,
      complete = T)
  print(p)
  dev.off()
}



