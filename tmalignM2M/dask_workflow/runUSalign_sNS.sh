#!/bin/bash
#To Run one to one pdb alignment using a semi-sequence independent method
#usage: sh runUSalign_sNS.sh <input_pdb> <run_against_pdb>

USALIGN_HOME=/path/to/tmalign/installed/bin/dir/USalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)

IN=$($USALIGN_HOME/USalign $inp_pdb $second_pdb -ter 1 -mm 6 -outfmt 2)
echo "$IN"

