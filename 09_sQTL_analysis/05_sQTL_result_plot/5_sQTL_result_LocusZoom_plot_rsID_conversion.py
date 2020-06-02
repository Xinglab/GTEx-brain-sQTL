#The purpose of this code is to map the SNPs found in our result to SNPs in the 1000 genome project based on chromosome + position
#SNPs that we used to calculate sQTL are from the GTEx genotype file and the name of those SNPs are something like this: 1_706368_A_G_b37
#However, in order to calculate LD between our SNPs and GWAS traits, we need to convert the sQTL SNPs into rs ID. This is what this code is doing. 
#The CEU.plink.*.tped files here contains positions of GWAS SNPs and their rs ID. We use this to convert our SNPs to rs ID (through genomic position)

import sys,glob
SNPs={}            #SNPs is a dictionary in which each key is the chromosome + position of a SNP and the value is empty
for l in open(sys.argv[1]):
	ls=l.strip().split('\t')[2].split('_')
	SNPs['\t'.join(ls[0:2])]=''      #different exons may have the same sQTL

matched={}       #the SNPs in our results that are also in the 1000 genome data
#fin_list=glob.glob(sys.argv[2]) # CEU.plink.\*.tped
fin_list=glob.glob('/path/to/tped/files/from/1000/genome/CEU.plink.'+sys.argv[2]+'.tped')


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
