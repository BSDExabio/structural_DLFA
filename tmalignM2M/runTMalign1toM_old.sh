#!/bin/bash
#To Run one to one pdb
#usage: sh runTmalign.sh <input_pdb> <run_against_pdb> <output_directory>
#To run one to many pdb
#usage: sh runTmalign.sh <input_pdb> <directory_of_pdb_to_run_against_many> <output_directory>
#The output will be saved in a file "output_directory/"
inp_pdb=$(readlink -f $1)
second_pdb=$(readlink -f $2)
outdir=$(readlink -f $3)

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

if [ ! -d $second_pdb ]; then
	IN=$(./TMalign $inp_pdb $second_pdb | grep 'RMSD\|TM-score' | head -2)
	rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) #$(echo $IN | cut -d "," -f 2)
        tmscore=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*') #$(echo $IN | cut -d "," -f -1)
        nm=$(basename $second_pdb)
	name=$(basename $inp_pdb)
	filename=$outdir/${name%.*}.txt
	header='PDB,RMSD,TM-score'
        echo $header > $filename
        echo "${nm%.*},$rmsd,$tmscore" >> $filename
	echo "PDB=${nm%.*}\tRMSD=$rmsd\tTM-score=$tmscore"

else
	header='PDB,RMSD,TM-score'
	name=$(basename $inp_pdb)
	#nm={$name%.pdb}
	filename=$outdir/${name%.*}.txt
	echo $header > $filename
	for pdb2 in $second_pdb/*.pdb; 
	do
		IN=$(./TMalign $inp_pdb $pdb2 | grep 'RMSD\|TM-score' | head -2)
		#echo $IN
		rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) #$(echo $IN | cut -d "," -f 2)
		tmscore=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*') #$(echo $IN | cut -d "," -f -1)
		#echo $rmsd
		#echo $tmscore
		nm=$(basename $pdb2)
		echo "${nm%.*},$rmsd,$tmscore" >> $filename

		#echo ${arrRMSD[1]}
	done

fi

echo "Output saved in $filename"
