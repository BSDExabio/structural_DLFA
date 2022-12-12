#!/bin/bash
#To Run one to one pdb alignment using a circular permutation-capable method
#Outputs log, rotation, and aligned structure files to the user-defined out_dir location
#usage:   sh runUSalign_cp_rbd.sh <input_pdb> <run_against_pdb> <output_directory_path>

USALIGN_HOME=/path/to/tmalign/installed/bin/dir/USalign

first_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)
out_dir=$(readlink -f $3)

$USALIGN_HOME/USalign $first_pdb $second_pdb -ter 1 -mm 3 -outfmt 0 -m $out_dir/trans_rot_matrix.dat > $out_dir/usalign.log	#-o $out_dir/USalign_cp 

