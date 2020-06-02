#make a sample file as this is required for rMATS
echo "SRR1473936.bam" > SRR1473936.sample.txt     #SRR1473936.bam is an example BAM file for one sample

#output folder for rMATS output
outdir="/output/folder"

#run rMATS
python /path/to/rMATS/rmats.py --b1 SRR1473936.sample.txt --od ${outdir} --tmp ${outdir}/.tmp \
  --anchorLength 1 --readLength 76 --gtf /path/to/GTF/GTF_file_name.gtf \
  -t paired --task prep --nthread 8 --statoff