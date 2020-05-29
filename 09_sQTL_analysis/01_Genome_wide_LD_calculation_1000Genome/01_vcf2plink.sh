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

input_file="1000_genome_genotype_by_chromosome.vcf.gz"
input_path="/path/to/1000/genome/genotype/data"
output_path="/output/path"
#export chrom="$SGE_TASK_ID"
export chrom="$2"
vcftools --gzvcf $input_file --plink-tped --out $output_path/$1.plink."$chrom"  --keep $input_path/$1.individuals.txt

#Input parameter:
#$1: Population (CEU/FIN/GBR/TSI/TRI)
#$2: chromosome number (1 to 22)

#Command parameter:
#--gzvcf: read compressed (gzipped) VCF files directly
#--plink-tped: outputs the genotype data in the PLINK transposed format with suffixes ".tped" and ".tfam"
#--keep: Provide files containing a list of individuals to include in subsequent analysis.
