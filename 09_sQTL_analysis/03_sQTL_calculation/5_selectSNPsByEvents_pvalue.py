import sys, numpy, glob
#variable1:full path + file name of association results of all exons
#e.g. Glimmps_each_exon_cis_Brain-Cortex/*
#variable2: full path + file name of all exon information
#e.g. /u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project/analysis/1_GLMM/sQTL/input_splicing/exon_info.fromGTF.SE.txt

#The purpose of this code is to select (significant) top SNP (with smallest p value. 
#if there are ties, we choose the one that is closest to the splice site)
# for each SE and output sQTL list

#example input:
#sys.argv=['nothing','Glimmps_each_exon_cis_Brain-FrontalCortexBA9/*','/u/nobackup/yxing/PROJECT/yidazhan/research/rotation_project/GTEx_brain_project_V7/analysis/6_sQTL_analysis/input_splicing/logit/JC/SE/exon_info.fromGTF.SE.txt'] 

fin_list= glob.glob(sys.argv[1])

def loadExonInfo(fin):             #generate a dictionary (key: exon ID. value: exon start position)
	exoninfo={}
	n=0
	for l in open(fin):
		if n==0:
			n+=1
			continue
		ls=l.strip().split('\t')
		exoninfo[ls[0]]=(int(ls[5])+int(ls[6]))/2 #middle of exon body (if we use python 2, this will be a integer, if we use python 3, this will be the exact number)
	return exoninfo

exoninfo=loadExonInfo(sys.argv[2])

fout=open('selected.SNPs.glmm.'+sys.argv[1].split('/')[0]+'.pvalue.txt','w')       #top SNPs for each exon
fout_sqtl=open('selected.sQTL.glmm.'+sys.argv[1].split('/')[0]+'.pvalue.txt','w')  #significant top SNPs (i.e., sQTL) for each exon

fin_list=glob.glob(sys.argv[1])
for fin in fin_list:          #for each exon
	#print fin
	se_id=fin.split('/')[1].split('.')[0]
	n=0
	events=[]       #store the p value of all SNPs
	events_snp=''   #store the event information of the top SNP
	events_pos=[]   #store the position of all SNPs
	events_pvalue=1 #store the p value of the top SNP
	events_dis=float("inf")   #store the distance of the top SNP
	for l in open(fin):       #for each SNP
		if n==0:
			n+=1
			continue     #skip the header
		ls=l.strip().split('\t')
		pvalue=float(ls[3])
		location=int(ls[2])
		dis=abs(exoninfo[se_id]-location)      #distance to the middle of the exon OR one of the two splice sites of the target exon (one splice site is the exon start and another splice site is the exons end)
		oridis=exoninfo[se_id]-location
		
		if events_snp=='':        #add the p value of the first SNP
			events_snp=l.strip()+'\t'+str(dis)+'\t'+str(oridis)
			events_pvalue=pvalue
			events_dis=dis
		else:          #for all the SNPs after the first one, we compare their p values with all the p values before them
			if pvalue<min(events):             #if we have a SNP with smaller p value
				events_snp=l.strip()+'\t'+str(dis)+'\t'+str(oridis)    #update the top SNP
				events_pvalue=pvalue
				events_dis=dis
			elif pvalue==min(events) and dis<events_dis:          #if we have ties on p value, we choose the one that is closer to the splice site
				events_snp=l.strip()+'\t'+str(dis)+'\t'+str(oridis)    #update the top SNP
				events_pvalue=pvalue
				events_dis=dis
		events.append(pvalue)
		events_pos.append(dis)
	fout.write(se_id+'\t'+events_snp+'\n')    #output the top SNP of the current exon
	if events_pvalue<10**-5:     #if the top SNP is also significant
		fout_sqtl.write(se_id+'\t'+events_snp+'\n')     #output the significant top SNP (i.e., sQTL) of the current exon
	
	
