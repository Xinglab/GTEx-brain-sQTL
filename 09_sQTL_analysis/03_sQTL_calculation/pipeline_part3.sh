#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=4G,h_rt=2:00:00
#$ -e ./part3_log/
#$ -o ./part3_log/

###This part of the code is to summarize the final result###

echo "pipeline_part3.sh started"
date

brain_region=$1
splicetype=$2
code_folder=$3
rootoutput=$4
counttype=$5
PSItype=$6

GWASdb="/path/to/GWAS_catalog/file/gwas_catalog_v1.0.1-associations_e89_r2017-07-31.tsv"
    
###########################
#8. summarize final result#
###########################
#based on raw p value
cd $rootoutput/$brain_region/"selected_sQTL_LD_GWAS_pvalue"
cat *.ld > "sig_SNPS_gwas.ld"
/usr/bin/python2.6 $code_folder/8_gwas_linkage_NHGRI.py $rootoutput/$brain_region/"selected_sQTL_LD_GWAS_pvalue"/"sig_SNPS_gwas.ld" $GWASdb

#based on permutation
cd $rootoutput/$brain_region/"selected_sQTL_LD_GWAS_permutation"
cat *.ld > "sig_SNPS_gwas.ld"
/usr/bin/python2.6 $code_folder/8_gwas_linkage_NHGRI.py $rootoutput/$brain_region/"selected_sQTL_LD_GWAS_permutation"/"sig_SNPS_gwas.ld" $GWASdb

echo "8_gwas_linkage_NHGRI.py finished"
date
