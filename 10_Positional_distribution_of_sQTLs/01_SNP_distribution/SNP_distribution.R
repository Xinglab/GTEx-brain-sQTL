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


#######################################################################################################
#Relationship between sQTL significance and distance of SNP to nearest splice site in 13 brain regions#
#######################################################################################################
print("part 5")
#read in the distance between sQTL and exon for each brain region
setwd(summaryinput)
distance=matrix(NA,nrow=0,5)
for (br in 1:length(brainregionlist)){
  brainregion=brainregionlist[br]
  formalbrainregion=formalbrainregionlist[br]
  
  output=try(suppressMessages(read.table(paste("all_info_exon_sQTL_GWAS_disease_in_",brainregion,"_",PSItype,"_",counttype,"_",splicetype,"_",type,".txt",sep=""),sep="\t",header=T)),silent=TRUE)   
  #each row is a exon-sQTL-GWAS trio, we need to change that to exon-sQTL pair
  
  if (!(inherits(output,"try-error"))){       #if there are significant sQTLs
    colname=c("exon_full_ID","genoposall","relative_distance","sQTL_p_value")
    significant.region=rep(formalbrainregion,dim(output)[1])
    temp=cbind(output[,colname],significant.region)
    distance=rbind(distance,unique(temp))
  }
}
colnames(distance)=c("exon","sQTL","distance","sQTL.p.value","brain.region")

library("plotrix")
library(ggsci)
library(scales)
library(ggplot2)
theme_ydz <- theme(#plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
  plot.title = element_text(size = rel(1.3), vjust = 2, hjust = 0.5, lineheight = 0.8),
  
  # axis
  axis.title.x = element_text(face="bold", size=16),
  axis.title.y = element_text(face="bold", size=16, angle=90, vjust=3),
  axis.text = element_text(size = rel(1.1)),
  axis.text.x = element_text(hjust = 0.5, vjust = 0, size=14),
  axis.text.y = element_text(vjust = 0.5, hjust = 0, size=14),
  axis.line = element_line(colour = "black"),
  
  # ticks
  axis.ticks = element_line(colour = 'black'),
  
  # legend
  legend.title=element_blank(),
  legend.text = element_text(size = rel(1.1),face="bold"),
  legend.position = "right",
  legend.direction = 'vertical',
  legend.background = element_blank(),
  legend.key = element_rect(fill = NA, colour = NA),
  
  # strip
  strip.text=element_text(size = rel(1.3)),
  
  # panel
  panel.background= element_blank(), 

  aspect.ratio=1,
  complete = T)

#plot#
setwd(outputpath)
filename=paste("distance_distribution_",PSItype,"_",counttype,"_",splicetype,"_",type,".pdf",sep="")
pdf(filename,width = 10, height = 6)
#brain region    
cbPalette <- c("#E53537",     
               "#EEA760",     
               "#C9E5C3",     
               "#9ACB3C",     
               "#028B45",     
               "#FFD923",     
               "#F0EC68",     
               "#6D67AF",     
               "#7F388D",    
               "#2999D5",     
               "#4EC3C7",     
               "#CD8ABC",     
               "#F287B7")     
COLOR="brain.region"    #this will control the legend name
temp=as.data.frame(distance)
brain.region=as.factor(as.character(as.matrix(distance[,"brain.region"])))
p=ggplot(temp,aes_string(x=(as.numeric(temp$distance)/1000), y=(-log10(as.numeric(temp$sQTL.p.value))), color=COLOR)) + 
  geom_point(size=2,alpha = 1) + 
  scale_color_manual(values=cbPalette) + 
  theme_ydz + 
  ggtitle(paste(splicetype,counttype,"\n")) +
  scale_x_continuous(limits = c(-200, 200)) + 
  scale_y_continuous(limits = c(0, ceiling(max(-log10(as.numeric(temp$sQTL.p.value)))/10)*10)) + 
  labs(x ="\nDistance to the nearest AS (kb)", y = expression(paste('-Log'[10], 'P-value')))
print(p)
dev.off()



