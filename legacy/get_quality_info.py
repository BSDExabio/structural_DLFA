import pandas as pd 
import argparse
from sadlsa_parsing_tools import SummaryParser as sp 
from sadlsa_parsing_tools import AlignmentParser as ap 
from sadlsa_parsing_tools import get_lengths_dic, extract_from_dic
import os 
parser = argparse.ArgumentParser()

parser.add_argument('-f', 
	dest='summary_file', 
	default='SummaryTable2.csv', 
	help='SAdLSA summary file')
parser.add_argument('-o', 
	dest='output_file', 
	default='sadlsa_alignment_quality.csv', 
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
parser.add_argument('--prob_cutoff', 
	dest='prob_cutoff', 
	default=0.8, 
	help='cutoff for high probability in 3A range for individual residues in alignment file')
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
enzonly = args.enzonly
cutoff = args.prob_cutoff



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



ldic = get_lengths_dic(lenlst)

master = []

for seq in seq_lst:
	aln_file = align_dir + '/' + seq + '/' + seq + '_aln.dat'
	start_lst, lines = ap.get_placements(aln_file)
	avg_qual = ap.get_avg_prob(start_lst, lines, num_align)
	prop_hi_3 = ap.get_prop_hiconf(lines, start_lst, num_align, cutoff, 3)
	aln_summary = ap.get_aln_summary(lines, start_lst, num_align, ldic, seq)
	seqlen = ldic[seq]
	for i in range(num_align):
		qual_lst = []
		qual_lst.append(seq)
		qual_lst.append(i+1)
		qual_lst.append(seqlen)
		qual_lst = extract_from_dic(aln_summary, qual_lst, i+1)
		qual_lst = extract_from_dic(avg_qual, qual_lst, i+1)
		qual_lst = extract_from_dic(prop_hi_3, qual_lst, i+1)
		master.append(qual_lst)


cols = ['SeqName', 'Alignment', 'SeqLength', 'tms1', 'tms2', 'SeqID', 'ProportionAligned', 'AvgProb3', 'AvgProb5', 'AvgProb9', 'PropHighProb']

df = pd.DataFrame(master, columns=cols)
df.to_csv(out_dir + '/' + output)






