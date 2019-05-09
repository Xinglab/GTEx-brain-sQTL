#!/bin/bash 
#$ -N plink2vcf
#$ -S /bin/bash                   
#$ -R y                                   
#$ -pe shared 1
#$ -l h_data=8G,h_rt=24:00:00
#$ -V       
#$ -cwd                                                                       
#$ -j y                                                                       
#$ -m be                                                  
#$ -M yidazhan@mail

#The purpose of this code is to change the file format from .vcf to .tped and .tfam for downstream calculation


#export chrom="$SGE_TASK_ID"
export chrom="$2"
vcftools --gzvcf ALL.chr"$chrom".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --plink-tped --out plink_allsnps/$1.plink."$chrom"  --keep /u/nobackup/yxing/NOBACKUP/sstein93/geuvadis/$1.individuals.txt

#Input parameter:
#$1: Population (CEU/FIN/GBR/TSI/TRI)
#$2: chromosome number (1 to 22)

#Command parameter:
#--gzvcf: read compressed (gzipped) VCF files directly
#--plink-tped: outputs the genotype data in the PLINK transposed format with suffixes ".tped" and ".tfam"
#--keep: Provide files containing a list of individuals to include in subsequent analysis.
