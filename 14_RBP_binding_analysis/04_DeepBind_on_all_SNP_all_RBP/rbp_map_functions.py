from pygr import seqdb
from collections import defaultdict
from pygr import seqdb
import math
import numpy as np
import cPickle as pickle
import matplotlib
#import matplotlib.pyplot as plt
from matplotlib import style

class exon(object):
	def __init__(self, line):
		(self.ID,
		 self.GeneID,
		 self.geneSymbol,
		 self.chrom,
		 self.strand,
		 self.exonStart_0base,
		 self.exonEnd,
		 self.upstreamES,
		 self.upstreamEE,
		 self.downstreamES,
		 self.downstreamEE,
		 self.ID,
		 self.IJC_SAMPLE_1,
		 self.SJC_SAMPLE_1,
		 self.IJC_SAMPLE_2,
		 self.SJC_SAMPLE_2,
		 self.IncFormLen,
		 self.SkipFormLen,
		 self.PValue,
		 self.FDR,
		 self.IncLevel1,
		 self.IncLevel2,
		 self.IncLevelDifference) = line.strip().split()
		self.line = line
		self.ID = int(self.ID)
		self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
		self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
		self.IJC_SAMPLE_2 = [int(x) for x in self.IJC_SAMPLE_2.split(",")]
		self.SJC_SAMPLE_2 = [int(x) for x in self.SJC_SAMPLE_2.split(",")]
		self.IncFormLen = int(self.IncFormLen)
		self.SkipFormLen = int(self.SkipFormLen)
		self.PValue = float(self.PValue)
		self.FDR = float(self.FDR)
		self.IncLevel1 =\
					  [float(x) for x in self.IncLevel1.split(",") if x != "NA"]
		self.IncLevel2 =\
					  [float(x) for x in self.IncLevel2.split(",") if x != "NA"]
		self.IncLevelDifference = float(self.IncLevelDifference)
		self.averagePsiSample1 = float(sum(self.IncLevel1)) / len(self.IncLevel2)
		self.averagePsiSample2 = float(sum(self.IncLevel2)) / len(self.IncLevel2)
		self.averageCountSample1 = float(sum(self.IJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1)) / float(len(self.IJC_SAMPLE_1))
		self.averageCountSample2 = float(sum(self.IJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2)) / float(len(self.IJC_SAMPLE_2))
		self.absDeltaPsi = math.fabs(self.IncLevelDifference)
		return
	
	def __str__(self):
		return self.line


def read_rMATS(fn):	
	#*--- filtering criteria ---*#
	# significanceFDR = 0.01  ## fdr siginificance level
	# significanceDeltaPSI = 0.1 ## abs(deltaPSI) significane level
	# backgroundExonsFDR = 0.5  ## fdr backgroundExons threshold
	# backgroundExonsCoverage = 50
	# foregroundExonsCoverage = 20
	# bgMinPSI = 0.15
	# bgMaxPSI = 0.85
	significanceFDR = 0.05  ## fdr siginificance level
	significanceDeltaPSI = 0.05 ## abs(deltaPSI) significane level
	backgroundExonsFDR = 0.5  ## fdr backgroundExons threshold
	backgroundExonsCoverage = 10
	foregroundExonsCoverage = 10
	bgMinPSI = 0.05
	bgMaxPSI = 0.95
	#*---- end filtering criteria ----*#
	up_exon_list = []
	dn_exon_list = []
	bg_exon_list = []
	with open(fn, 'r') as f:
		f.readline()
		for line in f:
			x = exon(line)
			if x.FDR<significanceFDR and x.averageCountSample1>foregroundExonsCoverage and x.averageCountSample2>foregroundExonsCoverage:
				if (x.IncLevelDifference > significanceDeltaPSI) :
					up_exon_list.append(x)
				elif (x.IncLevelDifference < -significanceDeltaPSI) :
					dn_exon_list.append(x)
			elif x.FDR>backgroundExonsFDR and \
				x.averageCountSample1 > backgroundExonsCoverage and\
                x.averageCountSample2 > backgroundExonsCoverage and\
                (x.averagePsiSample1>bgMinPSI and x.averagePsiSample1<bgMaxPSI) and\
                (x.averagePsiSample2>bgMinPSI and x.averagePsiSample2<bgMaxPSI) :
					bg_exon_list.append(x)
	exon_dict = {'Upregulated': up_exon_list, 'Downregulated': dn_exon_list, 'Background':bg_exon_list}
	return exon_dict


def read_motif_hexamers(fn, motif_len=6):
	motifs = set()
	with open(fn, 'r') as f:
		for line in f:
			line = line.strip()
			if not len(line)==motif_len: 
				raise Exception('Motif length not %s'%motif_len)
			motifs.add(line)
	return motifs


def read_genome(fn):
	'''the input of this function should be a FASTA file'''
	genome = seqdb.SequenceFileDB(fn)
	return genome


def rolling_window(seq, win_size):
	it = iter(seq)
	win = [it.next() for cnt in range(win_size)] # First window
	yield win
	for e in it: # Subsequent windows
		win[:-1] = win[1:]
		win[-1] = e
		yield win

def fetch_seq(genome, chr, start, end, strand):
	seq = ""
	try:
		seq = genome[chr][int(start):int(end)]
		if strand == "-":
			seq = -seq
	except:
		raise Exception('pygr cannot fetch sequences')
	return str(seq).upper()	


def compute_seq_motif_score(genome, chrom, start, end, strand, motif_set, motif_len=6, win_size=45):
	seq = fetch_seq(genome, chrom, start, end, strand)
	# occurences of 'any' of the hexamer sets
	found = [1 if ''.join(s) in motif_set else 0 for s in rolling_window(seq, motif_len)]
	motif_occurence = np.zeros(len(found))
	for i in range(len(found)):
		if found[i]==1:
			motif_occurence[i:(i+motif_len)] += 1
	motif_occurence = [1 if x>0 else 0 for x in motif_occurence]	
	
	# add trailing zero in the end to compensate for motif length
	#motif_occurence += [0]*(motif_len-1)
	
	# binning the occurences to get the window scores
	#motif_score = [1 if sum(w)>0 else 0 for w in rolling_window(motif_occurence, win_size)]
	#motif_score = [np.sum(w)/float(win_size) for w in rolling_window(motif_occurence, win_size)]
	motif_score = [np.sum(w)/float(win_size) for w in rolling_window(motif_occurence, win_size)]
	
	return motif_score

	

def RBP_map(exon_list, motif_set, genome, ext_out=250, ext_in=40, motif_len=6, win_size=45):
	#ext_out = 250
	#ext_in = 40
	start_padding = int((win_size-1)/2)
	end_padding = int((win_size-1)/2) + motif_len-1
	#end_padding = int((win_size-1)/2)
	upstream_score = np.zeros(ext_out+ext_in)
	dnstream_score = np.zeros(ext_out+ext_in)
	valid_seq_count = 0
	for x in exon_list:
		up_start, up_end = int(x.exonStart_0base) - ext_out - start_padding, int(x.exonStart_0base) + ext_in + end_padding
		dn_start, dn_end = int(x.exonEnd) - ext_in - start_padding, int(x.exonEnd) + ext_out + end_padding
		if x.strand == '-':
			up_start, dn_start = dn_start, up_start
			up_end, dn_end = dn_end, up_end
		
		this_up_score = compute_seq_motif_score(genome, x.chrom, up_start, up_end, x.strand, 
			motif_set, motif_len, win_size=45)
		this_dn_score = compute_seq_motif_score(genome, x.chrom, dn_start, dn_end, x.strand, 
			motif_set, motif_len, win_size=45)
		
		upstream_score += this_up_score
		dnstream_score += this_dn_score
		valid_seq_count +=1
	
	upstream_score /= float(valid_seq_count)
	dnstream_score /= float(valid_seq_count)
	return upstream_score, dnstream_score, valid_seq_count



def plot_map(score_dict, ext_out=250, ext_in=40, fn='rbp_map.pdf'):
	fig = plt.figure(figsize=(10,4))
	ax1 = fig.add_subplot(121) # for tgt-exon-upstream score
	ax1.set_title('TargetExon Upstream')
	ax2 = fig.add_subplot(122) # for tgt-exon-downstream score
	ax2.set_title('TargetExon Downstream')
	x_list = range(0, len(score_dict['Background'][0]))
	types = ['Upregulated', 'Downregulated', 'Background']
	colours = ['red', 'blue', 'black']
	ymax = 0
	for i in range(3):
		t = types[i]
		c = colours[i]
		ax1.plot(x_list, score_dict[t][0], c, label=t+', n='+str(score_dict[t][2]))
		ax2.plot(x_list, score_dict[t][1], c, label=t+', n='+str(score_dict[t][2]))
		this_max = np.max(np.concatenate([score_dict[t][0],score_dict[t][1]]))
		if this_max > ymax:
			ymax = this_max
	ymax = ymax * 1.2
	#print ymax
	ax1.set_ylim((0, ymax))
	ax2.set_ylim((0, ymax))
	ax1.set_xlim((0, ext_out+ext_in))
	ax2.set_xlim((0, ext_out+ext_in))
	ax1.axvline(ext_out, color='grey', linestyle='dashed')
	ax2.axvline(ext_in,  color='grey', linestyle='dashed')
	#ax1.legend()
	ax2.legend()
	ax1.set_xticklabels( ax1.get_xticks() - ext_out )
	ax2.set_xticklabels( ax2.get_xticks() - ext_in )
	plt.tight_layout()
	plt.savefig(fn)
	return

def return_seq_log_list(motif_info):
    motif_pos = []
    motif = motif_info.split('_')[-1]
    rbp = motif_info.split('_')[0]
    flag = False
    for x in  motif:
        if x == '[':
            motif_pos.append([])
            flag = True
            continue
        if x == ']':
            flag = False
            continue
        if flag:
            motif_pos[-1].append(x)
        else:
            motif_pos.append([x])
    all_pos_seq = []
    seq_list = generate_seq_list(motif_pos)
    return seq_list, rbp

def generate_seq_list(motif_pos):
    if len(motif_pos) == 1:
        return motif_pos[0]
    else:
        seq_list = []
        for motif in motif_pos[0]:
            for seq in generate_seq_list(motif_pos[1:]):
                seq_list.append(motif + seq)
        return seq_list

