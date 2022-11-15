#!/bin/bash
#To Run one to one pdb alignment using a circular permutation-capable method
#usage: sh runUSalign_cp.sh <input_pdb> <run_against_pdb>

### ISSUE: USalign's implimentation of Circular Permutation alignment now prints an extra line in the -outfmt 2 output
###        this results in needing a different parsing code than is currently implemented in the dask_tskmgr.py script
###        USalign has semi- and fully-sequence independent alignment methods; we should use those instead

USALIGN_HOME=/path/to/tmalign/installed/bin/dir/USalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)

IN=$($USALIGN_HOME/USalign $inp_pdb $second_pdb -ter 1 -mm 3 -outfmt 2)
echo "$IN"

