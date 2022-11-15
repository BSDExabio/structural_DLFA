#!/bin/bash
#To Run one to one pdb alignment using a semi-sequence independent method
#Outputs full log, rotation, and aligned structure files to the user-defined out_dir location
#usage:   sh runUSalign_sNS_save_aln.sh <input_pdb> <run_against_pdb> <output_directory_path>

USALIGN_HOME=/path/to/tmalign/installed/bin/dir/USalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)
out_dir=$(readlink -f $3)

$USALIGN_HOME/USalign $inp_pdb $second_pdb -ter 1 -mm 6 -outfmt 0 -o $out_dir/USalign_sNS -m $out_dir/trans_rot_matrix.dat > $out_dir/usalign.log

