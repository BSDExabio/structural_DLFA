#!/bin/bash
#To Run one to one pdb alignment using a fully-sequence independent method
#Outputs full log, rotation, and aligned structure files to the user-defined out_dir location
#usage:   sh runUSalign_fNS_save_aln.sh <input_pdb> <run_against_pdb> <output_directory_path>

USALIGN_HOME=/path/to/tmalign/installed/bin/dir/USalign

first_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)
out_dir=$(readlink -f $3)

$USALIGN_HOME/USalign $first_pdb $second_pdb -ter 1 -mm 5 -outfmt 0 -m $out_dir/trans_rot_matrix.dat -o $out_dir/USalign_fNS > $out_dir/usalign.log
rm $out_dir/*pml

