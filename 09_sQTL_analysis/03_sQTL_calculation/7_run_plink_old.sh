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
/u/local/apps/plink/1.08/plink --tfile /path/to/tped/files/from/1000/genome/CEU.plink.${chrom} --no-fid --no-parents --r2  --ld-snp-list $1 --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out $2/sig_SNPS_gwas_${chrom}


