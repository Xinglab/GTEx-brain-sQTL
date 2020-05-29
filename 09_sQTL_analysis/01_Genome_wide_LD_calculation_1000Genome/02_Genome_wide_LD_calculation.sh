#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=12G,h_rt=160:00:00,highp
#$ -e ./
#$ -o ./

#The purpose of this code is to calculate the LD of all pairwise SNPs (under the given filters) for each chromosome

outdir="/output/path"

chr=$1

#using 1000 Genome data:
plink --tfile /output/file/from/01_vcf2plink/by_chromosome --no-fid --no-parents --r2 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out $outdir/"LD_chr"$chr"_within_window"

#--tfile: the location of the .tped and .tfam file
#--no-fid: this allows you to use .fam or .ped files which lack family ID column
#--no-parents: this allows you to use .fam or .ped files which lack parental ID column

#LD statistic reports:
#--r2: calculates and reports squared correlations 
#--ld-window-kb: By default, when a limited window report is requested, every pair of variants with more than 1000 kilobases apart, is ignored. 
#--ld-window: By default, when a limited window report is requested, every pair of variants with at least 99999 variants between them, is ignored. 
#--ld-window-r2: With --r2, when a table format report is requested, pairs with r2 values less than 0.8 are filtered out of the report
