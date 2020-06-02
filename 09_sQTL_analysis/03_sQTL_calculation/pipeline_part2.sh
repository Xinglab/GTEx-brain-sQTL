#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=8G,h_rt=23:00:00
#$ -e ./part2_log/
#$ -o ./part2_log/

###This part of the code is to summarize glimmpse results and link the sQTLs with GWAS loci###

echo "pipeline_part2.sh started"
date

brain_region=$1
splicetype=$2
code_folder=$3
splicing_input=$4
rootoutput=$5
counttype=$6
PSItype=$7

cd $rootoutput  
#########################################################################
#5. select sQTL from glmm result (based on both p value and permutation)#
#########################################################################
cd $rootoutput/$brain_region
#based on raw p value
#example: python selectSNPsByEvents.py Glimmps_each_exon_cis/\* exon_info.fromGTF.SE.txt
/usr/bin/python2.6 $code_folder/5_selectSNPsByEvents_pvalue.py "Glimmps_each_exon_cis_"$brain_region"/*" $splicing_input/"exon_info.fromGTF."$splicetype".txt"
#based on permutation
source /u/local/Modules/default/init/modules.sh
module load R
Rscript $code_folder/5_selectSNPsByEvents_permutation.R "Glimmps_each_exon_cis_"$brain_region"/*" $splicing_input/"exon_info.fromGTF."$splicetype".txt"

echo "5_selectSNPsByEvents finished"
date

###############################################################################
#6. convert Genomic coordinates to rsID based on 1KGP CEU file (shayna folder)#
###############################################################################
cd $rootoutput/$brain_region
#based on raw p value
logfilepvalue="ID_conversion_for_"$brain_region".pvalue.log"
#example: python convertGPos2rsID.py ~/scratch/sQTL/selected.sQTL.glmm.Glimmps_each_exon_cis.txt
/usr/bin/python2.6 $code_folder/6_convertGPos2rsID.py $rootoutput/$brain_region/"selected.sQTL.glmm.Glimmps_each_exon_cis_"$brain_region".pvalue.txt" > $logfilepvalue
#based on permutation
logfilepermutation="ID_conversion_for_"$brain_region".permutation.log"
/usr/bin/python2.6 $code_folder/6_convertGPos2rsID.py $rootoutput/$brain_region/"selected.sQTL.glmm.Glimmps_each_exon_cis_"$brain_region".permutation.txt" > $logfilepermutation

echo "6_convertGPos2rsID.py finished"
date

#################
#7. calculate LD#
#################
cd $rootoutput/$brain_region
cp $code_folder/"chroms.txt" $rootoutput/$brain_region   #copy the chroms.txt file to the current brain region folder
mkdir "LD_log"     #this folder is used to store all the log files of LD calculation
#based on raw p value
jobnamepvalue="LD_"$brain_region"_pvalue"
inputSNPpvalue="selected.sQTL.glmm.Glimmps_each_exon_cis_"$brain_region".pvalue.txt.rsID.txt"
outdirpvalue="selected_sQTL_LD_GWAS_pvalue"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobnamepvalue -t 1-22:1 $code_folder/7_run_plink_old.sh $inputSNPpvalue $outdirpvalue
#based on permutation
jobnamepermutation="LD_"$brain_region"_permutation"
inputSNPpermutation="selected.sQTL.glmm.Glimmps_each_exon_cis_"$brain_region".permutation.txt.rsID.txt"
outdirpermutation="selected_sQTL_LD_GWAS_permutation"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobnamepermutation -t 1-22:1 $code_folder/7_run_plink_old.sh $inputSNPpermutation $outdirpermutation

echo "7_run_plink_old.sh finished"
date

###################################
#submit the third part of the code#
###################################
cd $code_folder
jobname="p3_"$PSItype"_"$counttype"_"$splicetype"_"$brain_region
#submit second part of the code when the second part is finished
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -hold_jid $jobnamepvalue,$jobnamepermutation -N $jobname $code_folder/pipeline_part3.sh $brain_region $splicetype $code_folder $rootoutput $counttype $PSItype

echo "pipeline_part3.sh submitted"
date
