#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=6G,h_rt=12:00:00
#$ -e ./part1_log/
#$ -o ./part1_log/

###This part of the code is to generate input and run glimmpse model to calculate correlation between exons and SNPs###

echo "pipeline_part1.sh started"
date

brain_region=$1
splicetype=$2
counttype=$6
PSItype=$7
code_folder=$3
splicing_input=$4
rootoutput=$5
#V7:
genotype_map_path="/path/to/genotype/file"
genotype_map_name="genotype_file_name_prefix"      #we don't use genotype information from X chromosome
genotype_vcf_path="/path/to/genotype/file"
genotype_vcf_name="name_of_genotype_vcf_file"
brain_metadata="/path/to/brain/tissue/annotation/file/gtex_v7_brain.txt"

mkdir -p $rootoutput

cd $rootoutput  
#################################################                 
#1. generate genotype file for each brain region#
#################################################
logfile="generate_genotype_file_for_"$brain_region".log"
#using V7 genotype
/usr/bin/python2.6 $code_folder/1_my_generateRaw_V7.py $splicing_input/$brain_region/"SRR_ID.txt" $genotype_map_path/$genotype_map_name $genotype_vcf_path/$genotype_vcf_name $brain_region $brain_metadata > $logfile

echo "1_my_generateRaw_V7.py finished"
date

###############################################
#2. split genotypes by chr for parallelization#
###############################################
cd $rootoutput/$brain_region
logfile="split_genotype_by_chr_for_"$brain_region".log"
/usr/bin/python2.6 $code_folder/2_split_raw_file_each_chr.py $rootoutput/$brain_region/$genotype_map_name"_tpose.raw" $genotype_map_path/$genotype_map_name $rootoutput/$brain_region > $logfile

echo "2_split_raw_file_each_chr.py finished"
date

#############################
#3. generate submission file#
#############################
cd $rootoutput/$brain_region
logfile="generate_submission_file_for_"$brain_region".log"
/usr/bin/python2.6 $code_folder/3_make.submission.py $splicing_input/"exon_info.fromGTF."$splicetype".txt" $splicing_input/$brain_region/"GTEx_brain_totalRC.txt" $splicing_input/$brain_region/"GTEx_brain_IC.txt" $brain_region $rootoutput/$brain_region $genotype_map_name $code_folder $PSItype $rootoutput > $logfile

#remove chrX as I did not include yet
cd $rootoutput/$brain_region/$brain_region
rm chrX -r
cd $rootoutput/$brain_region/$brain_region"_perm1"
rm chrX -r
cd $rootoutput/$brain_region/$brain_region"_perm2"
rm chrX -r
cd $rootoutput/$brain_region/$brain_region"_perm3"
rm chrX -r
cd $rootoutput/$brain_region/$brain_region"_perm4"
rm chrX -r
cd $rootoutput/$brain_region/$brain_region"_perm5"
rm chrX -r

echo "3_make.submission.py finished"
date

################
#4. submit jobs#
################
#create job.list for qsub job array (go to each permutation/original run folder)
cd $rootoutput/$brain_region/$brain_region
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname0="glimmps_"$brain_region
#example: qsub -t 1-43:1 /u/scratch2/p/panyang/sQTL/scripts/qsub.glimmps.sh
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname0 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh

cd $rootoutput/$brain_region/$brain_region"_perm1"
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname1="glimmps_"$brain_region"_perm1"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname1 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh

cd $rootoutput/$brain_region/$brain_region"_perm2"
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname2="glimmps_"$brain_region"_perm2"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname2 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh

cd $rootoutput/$brain_region/$brain_region"_perm3"
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname3="glimmps_"$brain_region"_perm3"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname3 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh

cd $rootoutput/$brain_region/$brain_region"_perm4"
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname4="glimmps_"$brain_region"_perm4"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname4 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh

cd $rootoutput/$brain_region/$brain_region"_perm5"
ls */*>job.list
jobnum=$(cat job.list | wc -l)
jobname5="glimmps_"$brain_region"_perm5"
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -N $jobname5 -t 1-$jobnum:1 $code_folder/4_qsub.glimmps.sh


echo "4_qsub.glimmps.sh finished"
date

####################################
#submit the second part of the code#
####################################
cd $code_folder
jobname="p2_"$PSItype"_"$counttype"_"$splicetype"_"$brain_region
#submit second part of the code when the first part is finished
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -hold_jid $jobname0,$jobname1,$jobname2,$jobname3,$jobname4,$jobname5 -N $jobname $code_folder/pipeline_part2.sh $brain_region $splicetype $code_folder $splicing_input $rootoutput $counttype $PSItype

echo "pipeline_part2.sh submitted"
date

