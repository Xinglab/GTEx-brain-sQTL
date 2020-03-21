#The difference between this code and Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_4_classification_and_plot.R
#is that in this code, instead of classifying exons purely based on the number of significant regions
#we use more stringent criteria:
#region specific: no more than 4 significant regions (p < 1e-5) & at least one region pass the FDR cutoff
#region ubiquitous: at least 10 regions pass the FDR cutoff


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

inputpath="/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot/Region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot_step_3_pick_SNP_and_collect_values_for_heatmap_more_stringent_version/single_job_run/result/200kb"

zscore=function(x){
  return((x-mean(x))/sd(x))     #z-score
  #return(x-mean(x))      #delta-psi
}

splicetype="SE"
type="pvalue"

source("/u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.3_motif_analysis/individual_exon_binding_peaks_scan/DeepBind/heatmap.3.R")

library(robust)
library(gplots)
library(ggplot2)
library(reshape)
library(RColorBrewer)
require(scales)
library(corrplot)
library(NMF)
library(fmsb)

pick_category=function(list_of_category,order_of_importance){      
  #given a list of categories and the order we want to use, this function will return one category from the list based on the order of importance
  reorderedlist=list_of_category[order(match(list_of_category, order_of_importance))]
  return(reorderedlist)
}

importance_order=c("dinucleotide", "SS", "exon", "<=300bp", ">300bp")

##################################################################################################
#read in the p value & beta matrix + number of significant region & location of SNP & rsID of SNP#
##################################################################################################
setwd(inputpath)
betamatrix200k=read.table("betamatrix_200kb.txt",sep="\t",header=T,check.names=F)
pvmatrix200k=read.table("pvmatrix_200kb.txt",sep="\t",header=T,check.names=F)
numsig.lc.rs200k=read.table("numsigregion_lc_rs_list200k.txt",sep="\t")


#########
#heatmap#
#########
#split the matrix into 2 groups
spergbeta=betamatrix200k[which(numsig.lc.rs200k[,"numsigregion.p"]<=4),]      #all the exons have at least one p value pass the FDR cutoff so we don't need to require that here
ubirgbeta=betamatrix200k[which(numsig.lc.rs200k[,"numsigregion.fdr"]>=10),]

colnames(spergbeta)=colnames(ubirgbeta)=formalbrainregionlist

#make a heatmap (z score transformed beta value)
for (group in 1:2){
  for (transformation in 1:2){
    if (group==1){
      #region specific group
      data2plot=spergbeta
      label1="Region_specific"
    }
    if (group==2){
      #region ubiquitous group
      data2plot=ubirgbeta
      label1="Region_ubiquitous"
    }
    
    if (transformation==1){
      #original beta value
      zbeta=data2plot
      label2="original_beta"
    }
    if (transformation==2){
      #z score transformation of the beta value
      zbeta=t(apply(data2plot,1,zscore))
      label2="z_score_transformed_beta"
    }
    
    setwd(inputpath)
    pdf(paste(label1,"_sQTL_beta_clustering_heatmap_",splicetype,"_",type,"_",label2,".pdf",sep=""),height=10,width=10)
    ownbreak = c(seq(range(zbeta)[1],-0.8,length=250),seq(-0.79,0.79,length=501),seq(0.8,range(zbeta)[2],length=250))
    heatmap.3(zbeta, 
              col = colorRampPalette(rev(brewer.pal(11,"RdBu")))(1000), 
              key = TRUE,
              breaks = ownbreak,
              Colv=FALSE)
    dev.off()
  }
}



#########################################
#statistical test of two group of events#
#########################################
specific=numsig.lc.rs200k[which(numsig.lc.rs200k[,"numsigregion.p"]<=4),]
ubiquitous=numsig.lc.rs200k[which(numsig.lc.rs200k[,"numsigregion.fdr"]>=10),]

#for the SNP location of each exon here, currently, the location is based on the location of the "top" SNP.
#the "top" SNP is also what we used for classification.
#however, for the statistical test here, we want to use the significant SNP with the best category (like what we did for SNP_to_SS_distance_vs_brain_region_specificity)

#read in the result from SNP_to_SS_distance_vs_brain_region_specificity
setwd("/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6.2_sQTL_SNP_annotation/result/SNP_to_SS_distance_vs_brain_region_specificity/closest_significant_SNP_categorized_distance/5_categories")
exonclass=read.table(paste(splicetype,"_",type,"_exonclass.txt",sep=""),sep="\t",header=T)[,1:length(brainregionlist)]

#update the location variable for the two classes of exons
newlclist=rep(NA,dim(specific)[1])
for (i in 1:dim(specific)[1]){
  currentexon=strsplit(rownames(specific)[i],split="~")[[1]][1]
  newlocation=pick_category(as.matrix(exonclass[currentexon,]),importance_order)[1]
  newlclist[i]=newlocation
}
specific=cbind(specific,newlclist)

newlclist=rep(NA,dim(ubiquitous)[1])
for (i in 1:dim(ubiquitous)[1]){
  currentexon=strsplit(rownames(ubiquitous)[i],split="~")[[1]][1]
  newlocation=pick_category(as.matrix(exonclass[currentexon,]),importance_order)[1]
  newlclist[i]=newlocation
}
ubiquitous=cbind(ubiquitous,newlclist)

#chi square test
data2test=cbind(table(specific[,"newlclist"]),table(ubiquitous[,"newlclist"]))
fisher.test(data2test,workspace=2e+07,hybrid=TRUE)
#For larger than  2 by 2 tables and 'hybrid = TRUE', asymptotic
#chi-squared probabilities are only used if the 'Cochran 
#conditions' are satisfied, that is if no cell has count zero, and 
#more than 80% of the cells have counts at least 5. 
#p-value < 2.2e-16 for (4,10),(3,11),(2,12) and (1,13)

#stacked bar plot
data2plot=as.data.frame(cbind(table(specific[,"newlclist"]),table(ubiquitous[,"newlclist"])))
colnames(data2plot)=c("specificpro","ubiquitouspro")
data2plot=melt(data2plot)
data2plot=cbind(data2plot,rep(names(table(specific[,"newlclist"])),2))
colnames(data2plot)=c("Group","Number.of.events","location")
data2plot$location <- factor(data2plot$location, levels = c("dinucleotide", "SS", "exon", "<=300bp", ">300bp"))
#c("Dinucleotide", "Splice Site", "Exon", "Intron (<=300bp)", "Intron (>300bp)")

setwd(inputpath)
pdf("stacked_bar_plot_proportion_comparison.pdf",height=4,width=4)
ggplot(data2plot, aes(x = Group, y = Number.of.events, fill = location)) + 
  geom_bar(position = "fill",width = 0.5,stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = scales::percent_format())+
  scale_fill_manual(values=c("#FDC086", "#BEAED4", "#0F0C73","#26908E","#626565"))+
  theme(
    # axis
    axis.title.x = element_text(face="bold", size=16, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(face="bold", size=12, margin = margin(t = 0, r = 10, b = 0, l = 0), angle=90),
    axis.text = element_text(size = rel(1.1)),
    axis.text.x = element_text(hjust = 0.5, vjust = 0, size=12),
    axis.text.y = element_text(vjust = 0.5, hjust = 0, size=12),
    axis.line = element_line(colour = "black"),
    
    #background
    panel.background = element_blank())
dev.off()



############
#radar plot#
############
radaroutputpath=paste(inputpath,"radar_plot",sep="/")
command=paste("mkdir -p",radaroutputpath)
system(command)

pair2plot=c("71949|ENSG00000073350.13_2|LLGL2|chr17|+|73570533|73570576|73570191|73570341|73570690|73570749~17_73561717_G_A_b37",
            "300088|ENSG00000135502.17_3|SLC26A10|chr12|+|58019018|58019113|58018905|58018939|58019192|58019247~12_58018979_C_T_b37",
            "88314|ENSG00000139793.18_3|MBNL2|chr13|+|98009049|98009103|97999294|97999321|98009735|98009889~13_97891187_T_C_b37")

pair2plot=c("9666|ENSG00000019144.18_3|PHLDB1|chr11|+|118484008|118484165|118478329|118478414|118484530|118484611~11_118478330_G_A_b37",
            "171698|ENSG00000205702.10_2|CYP2D7|chr22|-|42537634|42537685|42537271|42537349|42537877|42537934~22_42538103_T_C_b37")

pair2plot=c("359526|ENSG00000258366.7_3|RTEL1|chr20|+|62320854|62321001|62320407|62320485|62321102|62321218~20_62320968_T_C_b37",
            "263001|ENSG00000234127.8_3|TRIM26|chr6|-|30172432|30172542|30168811|30168915|30181081|30181142~6_30172513_A_C_b37",
            "258236|ENSG00000176731.11_2|C8orf59|chr8|-|86131464|86131592|86129614|86129731|86132534|86132613~8_86131463_A_T_b37",
            "258235|ENSG00000176731.11_2|C8orf59|chr8|-|86131464|86131589|86129614|86129731|86132534|86132582~8_86131463_A_T_b37",
            "270494|ENSG00000105136.20_3|ZNF419|chr19|+|58003480|58003579|57999332|58002965|58004223|58005151~19_58003580_A_G_b37",
            "222152|ENSG00000137312.14_2|FLOT1|chr6|-|30708966|30709110|30708455|30708575|30709390|30709481~6_30708955_C_T_b37",
            "344917|ENSG00000165915.13_2|SLC39A13|chr11|+|47434950|47435058|47433852|47434018|47435147|47435237~11_47434986_G_A_b37",
            "128125|ENSG00000171604.11_2|CXXC5|chr5|+|139059431|139059525|139028245|139028430|139059948|139060141~5_139059562_A_G_b37",
            "303254|ENSG00000184271.16_3|POU6F1|chr12|-|51593529|51593651|51592332|51592469|51597955|51598151~12_51599339_A_G_b37",
            "370678|ENSG00000077782.19_3|FGFR1|chr8|-|38286764|38286908|38285863|38285953|38287199|38287466~8_38286811_C_G_b37",
            "189436|ENSG00000186868.15_3|MAPT|chr17|+|44051750|44051837|44049224|44049311|44055740|44055806~17_43950441_C_A_b37")

#pair2plot=c("90646|ENSG00000122687.17_3|MRM2|chr7|-|2277649|2278572|2273929|2275199|2279052|2279342~7_2278226_C_T_b37",
##            "370678|ENSG00000077782.19_3|FGFR1|chr8|-|38286764|38286908|38285863|38285953|38287199|38287466~8_38286811_C_G_b37",
#            "359526|ENSG00000258366.7_3|RTEL1|chr20|+|62320854|62321001|62320407|62320485|62321102|62321218~20_62322699_A_G_b37",
#            "9912|ENSG00000169220.17_2|RGS14|chr5|+|176798475|176798590|176798350|176798397|176798873|176799167~5_176799992_T_C_b37",
#            "222152|ENSG00000137312.14_2|FLOT1|chr6|-|30708966|30709110|30708455|30708575|30709390|30709481~6_30708955_C_T_b37",
#            "229974|ENSG00000206530.10_3|CFAP44|chr3|-|113012798|113012885|113009699|113010595|113013533|113013668~3_113012797_G_A_b37",
#            "344917|ENSG00000165915.13_2|SLC39A13|chr11|+|47434950|47435058|47433852|47434018|47435147|47435237~11_47434986_G_A_b37",
#            "108171|ENSG00000035115.21_3|SH3YL1|chr2|-|224863|224920|218148|219001|229965|230044~2_221981_C_T_b37",
#            "108190|ENSG00000035115.21_3|SH3YL1|chr2|-|242797|242871|234159|234272|247537|247602~2_242793_G_C_b37",
#            "108192|ENSG00000035115.21_3|SH3YL1|chr2|-|243502|243562|234159|234272|247537|247602~2_242793_G_C_b37",
#            "128125|ENSG00000171604.11_2|CXXC5|chr5|+|139059431|139059525|139028245|139028430|139059948|139060141~5_139059562_A_G_b37",
#            "263001|ENSG00000234127.8_3|TRIM26|chr6|-|30172432|30172542|30168811|30168915|30181081|30181142~6_30172513_A_C_b37",
#            "303254|ENSG00000184271.16_3|POU6F1|chr12|-|51593529|51593651|51592332|51592469|51597955|51598151~12_51592761_C_A_b37",
#            "189436|ENSG00000186868.15_3|MAPT|chr17|+|44051750|44051837|44049224|44049311|44055740|44055806~17_44207360_C_T_b37")

for (pair in pair2plot){
  temp=strsplit(pair,split="~")[[1]]
  currentexon=temp[1]
  currentsnp=temp[2]
  genesymbol=strsplit(currentexon,split="\\|")[[1]][3]
  rsID=as.character(numsig.lc.rs200k[pair,"rslist"])
  
  data2plot=as.data.frame(matrix(0 , nrow=length(pair) , ncol=length(formalbrainregionlist)))
  rownames(data2plot)=pair
  colnames(data2plot)=formalbrainregionlist
  data2plot[1,]=pvmatrix200k[pair,]
  
  data2plot=-log10(data2plot)
  data2plot=rbind(rep(ceiling(max(data2plot)),length(formalbrainregionlist)) , rep(0,length(formalbrainregionlist)) , data2plot)
  
  numseg=ceiling(max(data2plot))     #number of segments for each axis (default 4). We want to show the significance cutoff so we need to change this for different events
  
  setwd(radaroutputpath)
  pdf(paste(pair,"_",currentsnp,"_",rsID,".pdf",sep=""),height=6,width=6)
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.7,0.5,0.1,0.9) )
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.7,0.5,0.1,0.4) )
  radarchart( data2plot  , axistype=1 , 
              title = paste(genesymbol,currentsnp,rsID,sep="\n"),
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(from=0,to=ceiling(max(data2plot)),length.out=numseg+1), cglwd=4,
              #custom labels
              vlcex=1,
              #other stuff
              centerzero=T,
              seg=numseg
  )
  #legend(x=1, y=1, legend = rownames(data2plot), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
  dev.off()
}



