#!/bin/bash

inp_pdb=$(readlink -f $1)
second_pdb=$(readlink -f $2) #either an pdb file or a directory of pdb files
outdir=$(readlink -f $3) #directory where the output will be stored
threshold=$(echo $4 | bc)

name=$(basename $inp_pdb)
#nm={$name%.pdb}
filename=$outdir/${name%.*}.txt

for pdb2 in $second_pdb/*; 
	do
        bsname=$(basename $pdb2)
        ext=${bsname##*.}
        echo "########################## $bsname extension is $ext %%%%%%%%%%%%%%%%%%%%"
        if [[ "$ext" != "pdb" || "$ext" != "cif" ]]; then
#        if [ "$ext" != "pdb" ]; then
            echo "$pdb2 EXTENSION FOUND IS A PDB FILE $ext"
	    #$ext=="pdb"
        fi
        if [ "$ext" != "pdb" ]; then 
		if  [ "$ext" != "cif" ]; then
			echo "EXTENSIONS dont match Skipping!!!"
		fi
	fi

        #if [ [ ! $ext=="pdb" ]  || [ ! $ext=="cif" ] ]; then
        #    echo "File extensions are not supported. Skipping!!!"
        #    continue
        #fi
done
