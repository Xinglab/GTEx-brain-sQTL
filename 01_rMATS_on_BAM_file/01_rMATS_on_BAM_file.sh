#make a sample file as this is required for rMATS
echo “SRR1473936.bam” > SRR1473936.sample.txt

#output folder for rMATS output
outdir=“/u/home/h/harryyan/bigdata/gtex”

#run rMATS
python /u/home/s/shiehshi/rMATS-2017-3-15/rmats.py --b1 SRR1473936.sample.txt --od ${outdir} --tmp ${outdir}/.tmp \
  --anchorLength 1 --readLength 76 --gtf /u/home/p/panyang/nobackup-yxing/references.annotations/37.chr/gencode.v26lift37.annotation.gtf \
  -t paired --task prep --nthread 8 --statoff