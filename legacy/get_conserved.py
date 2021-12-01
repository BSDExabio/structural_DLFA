import pandas as pd 
import argparse
from sadlsa_parsing_tools import SummaryParser as sp 
from sadlsa_parsing_tools import AlignmentParser as ap 
from sadlsa_parsing_tools import get_lengths_dic
import os 

parser = argparse.ArgumentParser()

parser.add_argument('-f', 
	dest='summary_file', 
	default='SummaryTable2.csv', 
	help='SAdLSA summary file')
parser.add_argument('-o', 
	dest='output_file', 
	default='sadlsa_alignment_conserved.csv', 
	help='name of output file containing alignment info')
parser.add_argument('--tms', 
	dest='tms_cutoff', 
	default='0.4', 
	help='TM score cutoff to choose sequences. Default=0.4')
parser.add_argument('-N', 
	dest='num_align', 
	default='5', 
	help='number of alignments to parse for each sequence. Default=5')
parser.add_argument('--wdir', 
	dest='working_dir', 
	default='NONE', 
	help='directory with summary files')
parser.add_argument('--alndir', 
	dest='align_dir', 
	default='NONE', 
	help='directory containing alignment files')
parser.add_argument('--outdir', 
	dest='output_directory', 
	default='NONE', 
	help='directory for output files')
parser.add_argument('--lenlst', 
	dest='len_lst', 
	default='de_hildenborough.lst', 
	help='file containing info on length for each sequence')
parser.add_argument('--cutoff', 
	dest='results_cutoff', 
	default=1, 
	help='Results cutoff - allows you to restrict results to just residues that appear at least N times (default = all)')
parser.add_argument('--enz', 
	dest='enzonly', 
	default='yes', 
	help='Allows you to get conserved residues for only sequences predicted to be enzymes. yes = enzymes only, no = all sequences. (default=yes)')


args = parser.parse_args()

file = args.summary_file
output = args.output_file
tms = float(args.tms_cutoff)
num_align = int(args.num_align)
working_dir = args.working_dir
align_dir = args.align_dir
out_dir = args.output_directory
lenlst = args.len_lst
cutoff = args.results_cutoff
enzonly = args.enzonly



if working_dir == 'NONE':
	working_dir = os.getcwd()

if align_dir == 'NONE':
	align_dir = working_dir + '/seq'

if out_dir == 'NONE':
	out_dir = working_dir + '/' + 'output'
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

if enzonly == 'yes':
	seq_lst = sp.get_enz_lst(file, tms)

if enzonly == 'no':
	seq_lst = sp.get_seq_lst(file, tms)


master_lst = []
for seq in seq_lst:
	conserved_lst = [[] for x in range(num_align + 1)]
	aln_file = align_dir + '/' + seq + '/' + seq + '_aln.dat'
	start_lst, lines = ap.get_placements(aln_file)
	conserved = ap.get_conserved(start_lst, lines, num_align)
	conserved = ap.how_conserved(conserved)
	for res in conserved:
		count = conserved[res]
		conserved_lst[count].append(res)
	conserved_lst[0].append(seq)
	master_lst.append(conserved_lst)

l = len(master_lst)


for r in range(l):
	for x in range(num_align+1):
		master_lst[r][x] = str(master_lst[r][x]).strip('[').strip(']')
		master_lst[r][x] = str(master_lst[r][x]).replace("'", "")

cols = []
cols.append('SeqName')
for i in range(num_align):
	cols.append(str(i+1) + ' conserved')

if not cutoff == 1:
	npop = int(cutoff) - 1
	for i in range(npop):
		for item in master_lst:
			del item[1]

if not cutoff == 1:
	npop = int(cutoff) - 1
	for i in range(npop):
		cols.pop(1)



df = pd.DataFrame(master_lst, columns=cols)
df.to_csv(output_dir + '/' + output)