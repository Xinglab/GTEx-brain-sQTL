outputpath="/path/to//12_region_ubiquitous_vs_specific_splicing_event_heatmap_radar_plot"
#####################
#check genotype swap#
#####################
setwd("/path/to/genotype/data")
command="cut -f 3,4,5 GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC_maf_0.01.vcf | sed 's/[\t]/\t/g' > Variant_ref_alt_lookup_table.txt"
system(command)

#R code start from here
VCFpath="/path/to/genotype/data"
Lookuptablename="Variant_ref_alt_lookup_table.txt"

setwd(VCFpath)
lookup=read.table(Lookuptablename,sep="\t",header=T)
rownames(lookup)=lookup[,"ID"]

get_ref_alt=function(x){
  snpid=as.matrix(x["ID"])
  ref=as.matrix(x["REF"])
  alt=as.matrix(x["ALT"])
  allele0=strsplit(snpid,split="_")[[1]][3]     #reference allele, the value is 0
  allele1=strsplit(snpid,split="_")[[1]][4]     #alternative allele, the value is 1
  if (ref==allele0 && alt==allele1){
    label="Correct"
  }else{
    label="Wrong"
  }
  return(label)
}

label=apply(lookup,1,get_ref_alt)

lookup=cbind(lookup,label)
setwd(outputpath)
write.table(lookup,"Genotype_swap_lookup_table.txt",sep="\t")




