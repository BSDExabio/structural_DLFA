#!/bin/bash
#To Run one to one pdb
#usage: sh runTMalign1toM.sh <input_pdb> <run_against_pdb> <output_directory> <tm-score_threshold>
#To run one to many pdb
#usage: sh runTMalign1toM.sh <input_pdb> <directory_of_pdb_to_run_against_many> <output_directory> <tm-score_threshold>
#The output will be saved in a file "output_directory/"

inp_pdb=$(readlink -f $1) 
second_pdb=$(readlink -f $2) #either an pdb file or a directory of pdb files
outdir=$(readlink -f $3) #directory where the output will be stored
threshold=$(echo $4 | bc)

echo "INPUT THRESHOLD=$threshold"
rerundir=$outdir/aligned_pdbs #if maximum TM-score is > than a threshold value then rerun TM-align and save the aligned PDBs here
rerunflag=0
header='PDB,RMSD,Length1,TM-score1,Length2,TM-score2,time'
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

if [ ! -d $second_pdb ]; then #if second PDB is a PDB file
	#Extract the RMSD and TM-score portions of the TMalign output
	IN=$(./TMalign $inp_pdb $second_pdb | grep 'RMSD\|TM-score\|CPU time' | head -8)
	#echo "****** $IN"
	rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) 
        tmscore1=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*' | head -1) 
        length1=$(echo $IN | grep -oP '(?<=LN=)[^ ]*' | tr -d , | head -1)
        tmscore2=$(echo $IN | grep -oP -m 2 '(?<=TM-score= )[^ ]*' | tail -1)
        length2=$(echo $IN | grep -oP -m 2 '(?<=LN=)[^ ]*' | tail -1 | tr -d ,)
        tim=$(echo $IN | grep -oP -m 2 '(?<=is )[^ ]*' | tail -1 )
        nm=$(basename $second_pdb)
	name=$(basename $inp_pdb)
	filename=$outdir/${name%.*}.txt
	#header='PDB,RMSD,TM-score'
        echo $header > $filename
        echo "${nm%.*},$rmsd,$length1,$tmscore1,$length2,$tmscore2" >> $filename
	echo "PDB=${nm%.*}\tRMSD=$rmsd\tLength1=$length1\tTM-score1=$tmscore1\tLength2=$length2\tTM-score2=$tmscore2\ttime=$tim"
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
	echo "MAXTMSCORE= $maxtmscore THRESHOLD=$threshold"
	rerunflag=$((`echo "$maxtmscore >= $threshold" | bc`))
	if [ $rerunflag -eq 1 ]; then
		echo "Since Maximum TM-score >= Threshold, rerunning the TMaling to save the aligned PDBs"
		if [ ! -d $rerundir ]; then
			mkdir -p $rerundir
		fi
		#Rerun TMalign and save the aligned pdb file 
		alignedfilename=${name%.*}_${nm%.*}
		if [ ! -d $rerundir/$alignedfilename ]; then
			mkdir -p $rerundir/$alignedfilename
		fi
		./TMalign $inp_pdb $second_pdb -o $rerundir/$alignedfilename/$alignedfilename\.pdb
		rerunflag=0
	fi
	

else #If second_pdb is a directory of PDBs
	name=$(basename $inp_pdb)
	#nm={$name%.pdb}
	filename=$outdir/${name%.*}.txt
	echo $header > $filename
	#Loop of all the PDB files in the directory of PDB files and run TMalign
	for pdb2 in $second_pdb/*.pdb; 
	do
		IN=$(./TMalign $inp_pdb $pdb2 | grep 'RMSD\|TM-score\|CPU time' | head -8)
		#echo "********************8##################### $IN"
                #echo "************************************ END *******************************"
		rmsd=$(echo $IN | grep -oP '(?<=RMSD= )[^ ]*' | tr -d ,) 
		tmscore1=$(echo $IN | grep -oP '(?<=TM-score= )[^ ]*' | head -1) 
		length1=$(echo $IN | grep -oP '(?<=LN=)[^ ]*' | tr -d , | head -1)
		tmscore2=$(echo $IN | grep -oP -m 2 '(?<=TM-score= )[^ ]*' | tail -1)
	        length2=$(echo $IN | grep -oP -m 2 '(?<=LN=)[^ ]*' | tail -1 | tr -d ,)
	        tim=$(echo $IN | grep -oP -m 2 '(?<=is )[^ ]*' | tail -1 )
                echo "******************* TIME = $tim"
		nm=$(basename $pdb2)
		echo "${nm%.*},$rmsd,$length1,$tmscore1,$length2,$tmscore2,$tim" >> $filename
		echo "PDB=${nm%.*}\tRMSD=$rmsd\tLength1=$length1\tTM-score1=$tmscore1\tLength2=$length2\tTM-score2=$tmscore2\ttime=$tim"
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
		echo "MAXTMSCORE= $maxtmscore THRESHOLD=$threshold"
		rerunflag=$((`echo "$maxtmscore >= $threshold" | bc`))
		
		if [ $rerunflag -eq 1 ]; then
			echo "Since Maximum TM-score >= $threshold Threshold, rerunning the TMaling to save the aligned PDBs"		
			if [ ! -d $rerundir ]; then
				mkdir -p $rerundir
			fi
			#Rerun TMalign and save the aligned pdb file 
			alignedfilename=${name%.*}_${nm%.*}
			if [ ! -d $rerundir/$alignedfilename ]; then
				mkdir -p $rerundir/$alignedfilename
			fi
			./TMalign $inp_pdb $pdb2 -o $rerundir/$alignedfilename/$alignedfilename\.pdb
			rerunflag=0

		fi
	done

fi

echo "Output saved in $filename"
if [ -d $rerundir ]; then
	echo "Aligned PDBs are saved in $rerundir"
fi
