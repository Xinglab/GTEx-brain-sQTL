#!/bin/sh

###
# tells the cluster to run the job from the current working directory
#$ -cwd

# tells the cluster what info to inform the user through email
#$ -m beas

# tells the cluster the user's email
#$ -M yidazhan@mail

# tells the cluster to use 1024MB of memory, and 1 hours of run time
#$ -l h_data=4G,h_rt=2:00:00

#$ -e ./11_1_sQTL_sashimi_plot_index_log/
#$ -o ./11_1_sQTL_sashimi_plot_index_log/
###


#step 1:
#get list of all BAM files
BAMfolder="/path/to/BAM/files"
cd $BAMfolder
BAMlist=(*.bam)

#step 2: for each BAM, we index it
for i in "${BAMlist[@]}"
do
	#echo "samtools index "$i
	echo "samtools index "$i > $i.sh    #generate the shell script
	qsub -V -cwd -m a -N $i -M yidazhan@mail -l h_data=6G,h_rt=6:00:00 $i.sh    #qsub the shell script
done


