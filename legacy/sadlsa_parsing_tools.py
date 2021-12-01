'''
Tools for parsing sadlsa outputs - looks for conserved residues in multiple alignments and gets alignment quality info
'''
from collections import Counter
import pandas as pd 

def get_lengths_dic(file): #intended for use on the de_hildenborough.lst file
	dic = {}
	f = open(file, 'r')
	lines = f.readlines()
	for line in lines:
		seq = line.split()[0]
		l = int(line.split()[1])
		dic[seq] = l 
	return dic 

def extract_from_dic(dic, lst, key):
	obj = dic[key]
	if isinstance(obj, list):
		for item in obj:
			lst.append(item)
	else:
		lst.append(obj)
	return lst



class AlignmentParser:

	def get_placements(file):
		'''
		Gets the locations of the start and finish of alignments in sadlsa alignment files.
		Also returns a list of lines from the alignment file
		'''
		start_lst = []
		f = open(file, 'r')
		lines = f.readlines()
		l = len(lines)
		for i in range(l):
			if lines[i][:3] == '###':
				start_lst.append(int(i))
		return start_lst, lines

	def get_conserved(start_lst, lines, num_alignments):
		'''
		Parses lines from alignment file and finds the conserved residues
		for top N alignments. Returns a dictionary with conserved residues for
		each alignment.
		'''
		dic = {}
		for i in range(num_alignments):
			lst = []
			for linenum in range(start_lst[i], start_lst[i+1]):
				line = lines[linenum].split()
				if '*' in line:
					lst.append(line[1])
			dic[i+1] = lst
		return dic

	def how_conserved(dic):
		'''
		takes dictionary from get_alignments, returns a dictionary with
		number of times a conserved residue is conserved
		'''
		lst = []
		for alignment in dic:
			conserved_lst = dic[alignment]
			for res in conserved_lst:
				lst.append(res)
		count = Counter(lst)
		return count  

	def get_avg_prob(start_lst, lines, num_alignments):
		dic = {}
		for i in range(num_alignments):
			lst3 = []
			lst5 = []
			lst9 = []
			for linenum in range(start_lst[i]+2, start_lst[i+1]-2):
				line = lines[linenum].split()
				lst3.append(float(line[-3]))
				lst5.append(float(line[-2]))
				lst9.append(float(line[-1]))
			mean3 = sum(lst3)/len(lst3)
			mean5 = sum(lst5)/len(lst5)
			mean9 = sum(lst9)/len(lst9)
			dic[i+1] = [mean3, mean5, mean9]
		return dic

	def get_prop_hiconf(lines, start_lst, num_alignments, cutoff, dist):
		d = {3:-3, 5:-2, 9:-1}
		col = d[dist]
		props = {}
		 
		for i in range(num_alignments):
			lst = []
			startline = start_lst[i]+2
			endline = start_lst[i+1]-2
			l = endline - startline
			for linenum in range(startline, endline):
				line = lines[linenum].split()
				prob = float(line[col])
				if prob >= cutoff:
					lst.append(prob)
			prop = len(lst)/l
			props[i+1] = prop  
		return props 

	def get_aln_summary(lines, start_lst, num_alignments, ldic, seqname):
		l = ldic[seqname]
		dic = {}
		for i in range(num_alignments):
			line = lines[start_lst[i]].split()
			tms1 = line[-3].strip('tms1=')
			tms2 = line[-2].strip('tms2=')
			sid = line[-1].strip('sid=').strip('%')
			naln = int(line[5].split('=')[1])
			prop = naln/l 
			dic[i+1] = [tms1, tms2, sid, prop]
		return dic 



class SummaryParser:

	def get_enz_lst(file, tms_cutoff):
	#get list of sadlsa genes predicted to be enzymes for parsing from summary table data
		lst = []
		df = pd.read_csv(file)
		l = df.shape[0]
		for i in range(l):
			if df.iloc[i]['DeepEC Binary'] == 'Enzyme':
				if df.iloc[i]['TM Score'] >= tms_cutoff:
					lst.append(df.iloc[i]['SeqName'])
		return lst

	def get_seq_lst(file, tms_cutoff):
		#get list of all genes for parsing from summary table
		lst = []
		df = pd.read_csv(file)
		l = df.shape[0]
		for i in range(l):
			if df.iloc[i]['TM Score'] >= tms_cutoff:
				lst.append(df.iloc[i]['SeqName'])

	def get_intersection(file, tms_cutoff, prob_cutoff):
		intersect = []
		sadlsa_good = []
		hhblits_good = []
		df = pd.read_csv(file)
		l = df.shape[0]
		for i in range(l):
			tms = df.iloc[i]['TM Score']
			prob = df.iloc[i]['HHblits Prob']
			seqname = df.iloc[i]['SeqName']
			if tms >= tms_cutoff and prob >= prob_cutoff:
				intersect.append([seqname, tms, prob])
			if tms >= tms_cutoff and prob < prob_cutoff:
				sadlsa_good.append([seqname, tms, prob])
			if tms < tms_cutoff and prob >= prob_cutoff:
				hhblits_good.append([seqname, tms, prob])
		return intersect, sadlsa_good, hhblits_good

	def get_matches(file):
		lst = []
		df = pd.read_csv(file)
		l = df.shape[0]
		for i in range(l):
			spdb = df.iloc[i]['SAdLSA PDB ID']
			spdb = spdb[:4]
			hhpdb = hhpdb[:4]
			hhpdb = df.iloc[i]['HHblits PDB ID']
			if spdb == hhpdb:
				seqname = df.iloc[i]['SeqName']
				lst.append([seqname, spdb])
		return lst



