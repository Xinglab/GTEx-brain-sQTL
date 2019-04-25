# variable1: full path + file name of significant sQTL file
# e.g. ~/scratch/sQTL/selected.sQTL.glmm.Glimmps_each_exon_cis.txt

#The purpose of this code is to map the SNPs found in our result to SNPs in the 1000 genome project based on chromosome + position
#SNPs that we used to calculate sQTL are from the GTEx genotype file and the name of those SNPs are something like this: 1_706368_A_G_b37
#However, in order to calculate LD between our SNPs and GWAS traits, we need to convert the sQTL SNPs into rs ID. This is what this code is doing. 
#The CEU.plink.*.tped files here contains positions of GWAS SNPs and their rs ID. We use this to convert our SNPs to rs ID (through genomic position)

#example input:
#sys.argv=['nothing','/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/sQTL_run/Brain-Cortex_copy/selected.sQTL.glmm.Glimmps_each_exon_cis_Brain-Cortex.txt','1']

import sys,glob
SNPs={}            #SNPs is a dictionary in which each key is the chromosome + position of a SNP and the value is empty
for l in open(sys.argv[1]):
	ls=l.strip().split('\t')[2].split('_')
	SNPs['\t'.join(ls[0:2])]=''      #different exons may have the same sQTL

matched={}       #the SNPs in our results that are also in the 1000 genome data
#fin_list=glob.glob(sys.argv[2]) # CEU.plink.\*.tped
fin_list=glob.glob('/u/nobackup/yxing/NOBACKUP/sstein93/1000GenomesVCF/VCF/plink_allsnps/CEU.plink.'+sys.argv[2]+'.tped')

#tped format:
#SNPs as rows and individuals as columns (The fifth and sixth fields are allele calls for the first sample in the .tfam file ('0' = no call); the 7th and 8th are allele calls for the second individual; and so on)
#The first four columns are from the MAP file (chromosome, SNP ID, genetic position, physical position), followed by the genotype data
#The reason why we change vcf file to tped is because it is faster to process. 
#The code that Shayna used to get the CEU.plink.*.tped  file is something like this:
#vcftools --gzvcf ALL.chr"$chrom".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --plink-tped --out plink_allsnps/$1.plink."$chrom"  --keep /u/nobackup/yxing/NOBACKUP/sstein93/geuvadis/$1.individuals.txt
#The samples used here is from CEU (European population). Using which sample will not influence anything here (here is just ID conversion) but will influence the LD calculation later.
#The vcf file from 1000 genome project of all samples from all populations can be downloaded here http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

fout=open(sys.argv[1]+'.rsID.txt','w')
fout1=open(sys.argv[1]+'.rsIDmap.txt','w')
for fin in fin_list:       #for each .tped file
	print fin
	for l in open(fin):     #for each SNP in the .tped file
		ls=l.strip().split('\t')
		if ls[0]+'\t'+ls[3] in SNPs:      #ls[0]+'\t'+ls[3] is the chromosome + position of a SNP. 
			matched[ls[0]+'\t'+ls[3]]=''    #if the SNP is also in our result, we add it to our results
			fout.write(ls[1]+'\n')      #output the rsID of the matched SNP
			fout1.write(ls[0]+'_'+ls[3]+'\t'+ls[1]+'\n')   #output the detailed information of the matched SNP (e.g., 3_60069\trs549251461\n)
for k in SNPs.keys():    #print out the SNPs that are not matched in 1000 genome
	if k not in matched:
		print k
