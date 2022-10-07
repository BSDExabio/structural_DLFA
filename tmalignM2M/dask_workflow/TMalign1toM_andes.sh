#!/bin/bash

# edit to point to your installation of TMalign
TMALIGN_HOME=/gpfs/alpine/proj-shared/bif135/rbd_work/TMalign

#To run one to many pdb
#usage: sh runTMalign1toM.sh <input_pdb> <directory_of_pdb_to_run_against_many> <output_directory> <tm-score_threshold>
#The output will be saved in a file "output_directory/"

inp_pdb=$(readlink -f $1) 	# path to a query pdb file
second_pdb=$(readlink -f $2) 	# path to a directory of library/target pdb files
outdir=$(readlink -f $3) 	# directory where the output will be stored
threshold=$(echo $4 | bc)	# floating point value between 0 and 1

header='PDB,RMSD,Length1,TM-score1,Length2,TM-score2,time,rerun'
rerundir=$outdir/aligned_pdbs #if maximum TM-score is > than a threshold value then rerun TM-align and save the aligned PDBs here
rerunflag=0

name=$(basename $inp_pdb)
filename=$outdir/${name%.*}.csv
echo $header > $filename
#Loop of all the PDB files in the directory of PDB files and run TMalign
for pdb2 in $second_pdb/*.pdb; 
do
	IN=$($TMALIGN_HOME/TMalign $inp_pdb $pdb2 | grep 'RMSD\|TM-score\|CPU time' | head -5)
	
	rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) 
	tmscore1=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*' | head -1) 
	length1=$(echo $IN | grep -oP '(?<=LN=)[^ ]*' | tr -d , | head -1)
	tmscore2=$(echo $IN | grep -oP -m 2 '(?<=TM-score= )[^ ]*' | tail -1)
        length2=$(echo $IN | grep -oP -m 2 '(?<=LN=)[^ ]*' | tail -1 | tr -d ,)
        tim=$(echo $IN | grep -oP -m 2 '(?<=is )[^ ]*' | tail -1 )

	nm=$(basename $pdb2)
	tm1=$(echo $tmscore1 | bc)
	tm2=$(echo $tmscore2 | bc)
	comp=$((`echo "$tm1 > $tm2" | bc`))
	compeq=$((`echo "$tm1 == $tm2" | bc`))
	if [ $compeq -eq 1 ]; then
		maxtmscore=$(echo $tmscore1 | bc)
	else
		if [ $comp -eq 1 ]; then
			maxtmscore=$(echo $tmscore1 | bc)
		else
			maxtmscore=$(echo $tmscore2 | bc)
		fi
			
	fi
	rerunflag=$((`echo "$maxtmscore >= $threshold" | bc`))

	echo "${nm%.*},$rmsd,$length1,$tmscore1,$length2,$tmscore2,$tim,$rerunflag" >> $filename
	
	if [ $rerunflag -eq 1 ]; then
		if [ ! -d $rerundir ]; then
			mkdir -p $rerundir
		fi
		#Rerun TMalign and save the aligned pdb file 
		alignedfilename=${name%.*}_${nm%.*}
		if [ ! -d $rerundir/$alignedfilename ]; then
			mkdir -p $rerundir/$alignedfilename
		fi
		$TMALIGN_HOME/TMalign $inp_pdb $pdb2 -o $rerundir/$alignedfilename/$alignedfilename > $rerundir/$alignedfilename/$alignedfilename.log
		rerunflag=0

	fi
done

if [ -d $rerundir ]; then
	tar -C $outdir -czf $outdir/aligned_pdbs.tar.gz aligned_pdbs/
	rm -Rf $rerundir
fi
