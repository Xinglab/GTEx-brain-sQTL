#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=2G,h_rt=1:00:00
#$ -e ./
#$ -o ./

code_folder="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/LD_calculation_using_GTEx_genotype_data"
outdir="/u/flashscratch/y/yidazhan/GTEx_V7_analysis/LD_calculation_using_GTEx_genotype_data/LD"

cd $outdir

for i in {1..22}       #each chromosome
do 
jobname="chr_"$i"_LD"
#echo $jobname
/u/systems/UGE8.6.4/bin/lx-amd64/qsub -N $jobname $code_folder/0_LD_calculation.sh $i
done

