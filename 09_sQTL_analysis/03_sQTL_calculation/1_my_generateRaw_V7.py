import sys,numpy,os
#What the code is doing:
#1. load all the sample information + count the number of donors + pair the SRR ID with donor ID
#2. load the samples in the current brain region
#3. load the .map file (information about all the SNPs after the filter)
#4. parse the vcf file. This is the key purpose of this code. The idea is that the vcf file contains information for all samples and all SNPs. 
#But we are just running our analysis in each brain region separately. So we don't want to load all the information every time. 
#So we just want to get the genotype (SNPs after filter) of samples in the current brain region (and also in the vcf file because some samples may not have genotype information)


def loadSampleDict(fin):     #this function pairs SRR ID with donor ID and count the number of donors in the SRR run table
	sample_dict={}
	donor_counts={}
	n=0
	for l in open(fin):
		if n==0:
			n+=1
			continue
		ls=l.strip().split('\t')
		sample_dict[ls[15]]='-'.join(ls[17].split('-')[0:2])     #my modification
		donor_counts['-'.join(ls[17].split('-')[:2])]=''         #my modification
	print 'total donors:', len(donor_counts)
	return sample_dict

def loadSampleList(fin):          #this function returns the SRR ID of samples in the current brain region
	sample_list=[]
	for l in open(fin):
			sample_list.append(l.strip())
	return sample_list

def loadMAPfile(fin):         #this function loads the .map file into a dictionary: each key is the SNP ID and the value is ' '
	snp_dict={}
	for l in open(fin):
		ls=l.strip().split('\t')
		snp_dict[ls[1]]=''
	print 'total snp counts:',len(snp_dict)
	return snp_dict

def parseVCFfile(fin, filename, genotypename, snp_dict, sample_list, sample_dict):
	genotype_convert={'0/0':'0','1/1':'2','1/0':'1','0/1':'1','./.':'NA'}      #convert genotype into number
	header=[]   #all the donor ID
	sample_list_has=[]       #samples that have genotype information
	os.system('mkdir -p '+genotypename)
	fout=open(genotypename+'/'+filename.split('/')[-1]+'_tpose.raw','w')
	for l in open(fin):          #for each line (the vcf file can be grouped into 3 parts)  
		if l.startswith('##'):     #part 1:file info part, we do nothing.
			continue
		if l.startswith('#CHROM'):    #part 2: the #CHROM line, we get donor ID and phenotype information.
			ls=l.strip().split('\t')
			header=map(lambda x: x.split('_')[0],ls[9:])   #get shorter donor ID     (different from V6)
			row=dict(zip(header, ls[9:]))      #pair shorter donor ID with full donor ID: 'GTEX-P4PP':GTEX-P4PP-0004-SM-2H2WW
			line, line2, line34, line5, line6=['']*5      #these are used to store the phenotype information. However, we don't really have phenotype information here so we just give them arbitrary values
			for s in sample_list:
				if sample_dict[s] in row:   #only samples with genotype information
					line+='\t'+s
					line2+='\t1'
					line34+='\t0'
					line5+='\t1'
					line6+='\t-9'
					sample_list_has.append(s)
				else:      #if the sample doesn't have genotype information
					print s, sample_dict[s]
			fout.write('FID'+line+'\nIID'+line2+'\nPAT'+line34+'\nMAT'+line34+'\nSEX'+line5+'\nPHENOTYPE'+line6+'\n')    #output all the phenotype information of samples that have genotype information
			continue
		ls=l.strip().split('\t')      #part 3: the SNP part, we get the SNP information of samples with genotype.
		if ls[2] not in snp_dict:     #there may be SNPs that in the vcf file but not in the .map file because we filtered out some SNPs when we generate the .map file
			continue
		row=dict(zip(header, ls[9:]))   #pair shorter donor ID with genotype information, e.g., 'GTEX-145MO': '0/1'
		line=ls[2]     #SNP ID
		for s in sample_list_has:    #get the genotype information for all the samples that have genotype information
			line+='\t'+genotype_convert[row[sample_dict[s]]]        #            (different from V6)
		fout.write(line+'\n')              #output all the genotype information of the current SNP
	print 'sample with genotype:',len(sample_list_has)

def main():

	brain_metadata=sys.argv[5]
	sample_dict=loadSampleDict(brain_metadata) # Sra run table  #my modification
	print 'sample and donor list:', len(sample_dict)
	sample_list=loadSampleList(sys.argv[1]) 
	print 'sample list:',len(sample_list)
	snp_dict=loadMAPfile(sys.argv[2]) #map file
	parseVCFfile(sys.argv[3], sys.argv[2], sys.argv[4], snp_dict, sample_list, sample_dict) #vcf file, map file,...
	print 'done.'

if __name__ == '__main__':
	main()


