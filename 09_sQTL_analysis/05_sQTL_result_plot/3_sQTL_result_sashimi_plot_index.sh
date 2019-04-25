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



###pay attention: index all the BAM files first before running sashimi plot. This only needs to be done once for each BAM file###
###the command for indexing BAM file is "samtools index bam.bam", the output index file will be in the same location of the BAM files###
###for incomplete BAM file, we can use "samtools view SRRXXXXXXX.bam | tail | hexdump -C" to double check
###see http://seqanswers.com/forums/showthread.php?t=15363

#step 1:
#get list of all BAM files
BAMfolder="/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/data/V7_all_bam_soft_link"
cd $BAMfolder
BAMlist=(*.bam)

#step 2: for each BAM, we index it
for i in "${BAMlist[@]}"
do
	#echo "samtools index "$i
	echo "/u/nobackup/yxing/PROJECT/yidazhan/research/software/samtools-1.5/bin/samtools index "$i > $i.sh    #generate the shell script
	/u/systems/UGE8.0.1/bin/lx-amd64/qsub -V -cwd -m a -N $i -M yidazhan@mail -l h_data=6G,h_rt=6:00:00 $i.sh    #qsub the shell script
done


