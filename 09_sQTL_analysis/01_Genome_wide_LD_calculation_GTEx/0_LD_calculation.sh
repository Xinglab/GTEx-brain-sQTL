#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=12G,h_rt=160:00:00,highp
#$ -e ./
#$ -o ./

#The purpose of this code is to calculate the LD of all pairwise SNPs (under the given filters) for each chromosome
#so that later if we want to know the r square between any two SNPs on a given chromosome, we can directly check the result here without doing any calculation
#This code is independent of the sQTL pipeline
#The result is used in the haploview plot part


outdir="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/LD_calculation_using_GTEx_genotype_data/LD"

chr=$1

#using GTEx data (we haven't run this yet so we don't have this result):
/u/local/apps/plink/1.08/plink --tfile /u/project/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/data/raw_data/files/genotype/V7_whole_exon_sequencing_modified/split_chromosome/vcf2rsID_dbSNP147_GRCh37p13.chr${chr} --no-fid --no-parents --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out $outdir/"LD_chr"$chr"_within_window"


