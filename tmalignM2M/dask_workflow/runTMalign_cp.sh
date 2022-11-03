#!/bin/bash
#To Run one to one pdb
#usage: sh runTMalign_cp_rbd.sh <input_pdb> <run_against_pdb>

TMALIGN_HOME=/path/to/tmalign/installed/bin/dir/TMalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)

IN=$($TMALIGN_HOME/TMalign $inp_pdb $second_pdb -cp -ter 1 -outfmt 2)
echo "$IN"

