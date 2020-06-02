#!/bin/sh

#$ -cwd 
#$ -m beas
#$ -M yidazhan@mail
#$ -l h_data=2G,h_rt=1:00:00
#$ -e ./
#$ -o ./

code_folder="/path/of/code/02_Genome_wide_LD_calculation.sh"
outdir="/output/path"

cd $outdir

for i in {1..22}       #each chromosome
do 
jobname="chr_"$i"_LD"
#echo $jobname
qsub -N $jobname $code_folder/02_Genome_wide_LD_calculation.sh $i
done

