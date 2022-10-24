#!/bin/bash
#To Run one to one pdb
#usage: sh runTMalign1to1_cp.sh <input_pdb> <run_against_pdb>

TMALIGN_HOME=/gpfs/alpine/proj-shared/bif135/rbd_work/TMalign

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2)

#Extract the RMSD and TM-score portions of the TMalign output
IN=$($TMALIGN_HOME/TMalign $inp_pdb $second_pdb -cp | grep 'RMSD\|TM-score\|CPU time' | head -8)
nAln=$(echo $IN | grep -oP '(?<=Aligned length= )[^ ]*' | tr -d ,) 
rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) 
tmscore1=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*' | head -1) 
length1=$(echo $IN | grep -oP '(?<=LN=)[^ ]*' | tr -d , | head -1)
tmscore2=$(echo $IN | grep -oP -m 2 '(?<=TM-score= )[^ ]*' | tail -1)
length2=$(echo $IN | grep -oP -m 2 '(?<=LN=)[^ ]*' | tail -1 | tr -d ,)
tim=$(echo $IN | grep -oP -m 2 '(?<=Total CPU time is )[^ ]*')
#tim=$(echo $IN | grep -oP -m 2 '(?<=is )[^ ]*' | tail -1 )
echo "$rmsd $nAln $length1 $tmscore1 $length2 $tmscore2 $tim"

