#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=2G,h_rt=1:00:00
#$ -e ./
#$ -o ./

code_folder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/scripts"
outdir="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/LD"

cd $outdir

for i in {1..22}       #each chromosome
do 
jobname="chr_"$i"_LD"
#echo $jobname
/u/systems/UGE8.0.1/bin/lx-amd64/qsub -N $jobname $code_folder/Genome_wide_LD_calculation.sh $i
done

