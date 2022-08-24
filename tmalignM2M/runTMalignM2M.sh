#!/bin/bash
#Usage: sh ./runTMalignM2M.sh <input_dir_of_pdbs> <directory_of_against_pdb_folder> <output_dir> <tmscore_cutoff_threshold>

inp_dir=$(readlink -f $1)
against_pdb_folder=$(readlink -f $2)
outdir=$(readlink -f $3)
threshold=$4

root_dir=$(dirname $(readlink -f $0))

#TODO: Add validation for parameters

if [ ! -d $inp_dir ]; then
    echo "$inp_dir is not a directory"
    exit 1
fi

if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

for pdb in $inp_dir/*.pdb; 
do
	pdbname=$(basename $pdb)
	outname=$outdir/${pdbname%.pdb}
	sh $root_dir/runTMalign1toM.sh $pdb $against_pdb_folder $outname $threshold
done



