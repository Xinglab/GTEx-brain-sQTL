#!/bin/bash 
#$ -S /bin/bash                   
#$ -R y                                   
#$ -l h_data=6G,h_rt=3:00:00
#$ -V       
#$ -cwd                                                                       
#$ -j y                                                                       
#$ -m a                                                  
#$ -M yidazhan@mail  
#$ -e ./LD_log/
#$ -o ./LD_log/    
# 


# chroms.txt is list of chromosomes, one chromosome per line (1-22,X) 
export chrom=`sed -n ${SGE_TASK_ID}p chroms.txt`
echo $chorm
echo "input SNP list", $1
echo "out dir", $2
mkdir -p $2
/u/local/apps/plink/1.08/plink --tfile /u/nobackup/yxing/NOBACKUP/sstein93/1000GenomesVCF/VCF/plink_allsnps/CEU.plink.${chrom} --no-fid --no-parents --r2  --ld-snp-list $1 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out $2/sig_SNPS_gwas_${chrom}


#--tfile {prefix} : Specify .tped + .tfam filename prefix (default 'plink').
#--no-fid         : .fam/.ped file does not contain column 1 (family ID).
#--no-parents     : .fam/.ped file does not contain columns 3-4 (parents).
#--r2       : LD statistic reports.  --r yields raw inter-variant correlations, while --r2 reports their squares.
#--ld-snp-list [f]  : Restrict first --r/--r2 var. to those named in the file.
#--ld-window-kb [x] : Set --r/--r2 max kb pairwise distance (usually 1000).
#--ld-window [ct+1] : Set --r/--r2 max variant ct pairwise distance (usu. 10).
#--ld-window-r2 [x] : Set threshold for --r2 report inclusion (usually 0.2).
#--out [prefix]   : Specify prefix for output files.