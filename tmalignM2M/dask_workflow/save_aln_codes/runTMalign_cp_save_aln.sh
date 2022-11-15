#!/bin/bash
#To Run one to one pdb alignment using a circular permutation-capable method
#Outputs log, rotation, and aligned structure files to the user-defined out_dir location
#usage:   sh runTMalign_cp_rbd.sh <input_pdb> <run_against_pdb> <output_directory_path>

TMALIGN_HOME=/path/to/tmalign/installed/bin/dir/TMalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)
out_dir=$(readlink -f $3)

$TMALIGN_HOME/TMalign $inp_pdb $second_pdb -cp -ter 1 -outfmt 0 -o $out_dir/TMalign_cp -m $out_dir/trans_rot_matrix.dat > $out_dir/tmalign.log

