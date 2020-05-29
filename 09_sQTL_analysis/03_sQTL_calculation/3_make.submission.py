import sys, os
#the purpose of this code is to devide exons into batches and generate the shell script for each batch (this is under each permutation/original run and also each chromosome in each permutation/original run)

batch=500       #500 exons a time

out_suffix= sys.argv[4] # brain region 
exon_info_fin=sys.argv[1]
DP_fin=sys.argv[2] #GTEx_brain_totalRC.txt
IC_fin=sys.argv[3] #GTEx_brain_IC.txt
genotype_name= sys.argv[5] #genotype folder
map_prefix=sys.argv[6].split('.map')[0]         #e.g. Genotype_450_chr122_maf0.05
code_folder=sys.argv[7]
PSItype=sys.argv[8]
rootoutput=sys.argv[9]
task_list=['no_perm','perm1','perm2','perm3','perm4','perm5']    #original run + 5 permutations

allchr=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']


for perm_n in task_list:     #for each round
	chr_dist={}        #count the number of exons in each chromosome
	n=0
	for l in open(exon_info_fin): #exon_info.fromGTF.SE.txt 
		if n==0:
			n+=1
			continue       #skip the first row
		ls=l.strip().split('\t')
		if ls[3] in allchr:       #if the chromosome is 1-22 or X (exons from chrY and other chromosomes are ignored)
			if ls[3] in chr_dist:
				chr_dist[ls[3]]+=1            #count the number of exons on each chromosome
			else:
				chr_dist[ls[3]]=1
	for k in chr_dist.keys():         #for each chromosome
		print k, chr_dist[k]
		os.system('mkdir -p '+out_suffix)     #generate the folder for the original run
		
		if perm_n!='no_perm':
			os.system('mkdir -p '+out_suffix+'_'+perm_n+'/'+k)     #generate the folder for each permutation and also generate folders for each chromosome under the permutation folder
		else:
			os.system('mkdir -p '+out_suffix+'/'+k)      #generate folders for each chromosome under the original run folder
		n=0
		n=chr_dist[k]/batch      #n=the number of batches for the k chromosome (if n is multiple of the batch size); chr_dist[k] is the number of exons in the k chromosome
		
		if chr_dist[k]%batch!=0:     #if n is not multiple of the batch size, n=n+1 (now n=number of batches)
			n+=1
			
		
		for i in range(0,n):   #for each batch, we generate the shell script
			if perm_n=='no_perm':
				fout=open(out_suffix+'/'+k+'/'+k+'.glimmps.'+str(i)+'.sh','w')
				fout.write('source /u/local/Modules/default/init/modules.sh'+'\nmodule load R'+'\ncd '+rootoutput+'/'+out_suffix+'/\nRscript '+ code_folder + '/sQTLregress.oneexon.cis.R '+exon_info_fin+' '+DP_fin+' '+IC_fin+' '+map_prefix+' '+k+' '+str(i*batch+1)+' '+str((i+1)*batch)+' '+out_suffix+' '+genotype_name+' '+PSItype+'\n')   #my code
				fout.close()
			else:
				fout=open(out_suffix+'_'+perm_n+'/'+k+'/'+k+'.glimmps.'+perm_n+'.'+str(i)+'.sh','w')
	                        fout.write('source /u/local/Modules/default/init/modules.sh'+'\nmodule load R'+'\ncd '+rootoutput+'/'+out_suffix+'/\nRscript ' + code_folder + '/sQTLregress.oneexon.cis.perm.R '+exon_info_fin+' '+DP_fin+' '+IC_fin+' '+map_prefix+' '+k+' '+str(i*batch+1)+' '+str((i+1)*batch)+' '+perm_n.split('perm')[-1]+' '+out_suffix+' '+genotype_name+' '+PSItype+'\n')   #my code
	                        fout.close()
