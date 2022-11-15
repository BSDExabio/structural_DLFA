#!/bin/bash
inpdir=$(readlink -f $1)
truepdb_dir=$(readlink -f $2)
outdir=$(readlink -f $3)
if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

filename=$outdir/tmalignout.txt

for file in $inpdir/*.pdb;
do
	filenm=$(basename $file)
	name=${filenm%.*}
	IN=$(./TMalign $file $truepdb_dir/$name* | grep 'RMSD\|TM-score')
	echo $IN
done
