(These are all files original in the top-level directory that have been 
moved out of the way since they're largely deprecated.)

Scripts to get alignment information including alignment quality and conserved 
residues for SAdLSA outputs.

* `extract_aln.sh` - script that extracts and untars sadlsa alignments. Expects 
  one argument, the name of the directory where the directories for each 
  sequence's sadlsa alignments are located. Right now extracts only the pdb70 alignments. 

* `get_conserved.py` - script that gets conserved residues for the top N number 
  of alignments and outputs to a csv file. See output folder for the output from 
  this script. Default outputs results for hypothetical proteins predicted 
  to be enzmyes with a TM score 1 of >= 0.4.  Many parameters are customizable, see 
  script for all available options and defaults. Requires the (What?)

   > Ex: `python get_conserved.py -f SummaryTable2021.csv -o conserved_residues.
 csv --tms 0.4 -N 5`

* `get_quality_info.py` - Extracts the info of the quality of the alignment 
   from sadlsa alignments. Outputs a csv with the TM score 1, TM score 2, 
   sequence identity to the target, proportion aligned, the average 
   probability of that two aligned residues are within X angstroms (3, 5, and 9 
   A distance bins), and the proportion of residues with a high probability that 
   the aligned residues are within 3 A for the top N alignments for each sequence. 
   There are also many customizable parameters in this scripts. Default outputs results for 
   proteins predicted to be enzymes with a tms1 >= 0.4. See the output folder for an 
   output from this script.

   > Ex: `get_quality_info.py -f SummaryTable2021.csv -o alignment_quality_info.
csv -N 5 --tms 0.4`
