#!/bin/bash

while getopts d: flag
do
	case "${flag}" in 
		d) directory=${OPTARG};;
	esac
done

echo $directory

if [[ ! -e seq ]]; then
	mkdir seq
fi

d=$PWD

for f in $directory/*
do
	seq=${f##*/}


	aln=${directory}/${seq}/sadlsa_pdb70_210310/${seq}_aln.dat.gz
	aln_unzipped=${directory}/${seq}/sadlsa_pdb70_210310/${seq}_aln.dat

	if [[ -f "$aln" ]]; then
		
		
		if [[ ! -e seq/${seq} ]]; then
			mkdir seq/$seq
		fi

		cp $aln ${PWD}/seq/$seq
		gzip -d ${PWD}/seq/$seq/${seq}_aln.dat.gz
	fi


	if [[ -f "$aln_unzipped" ]]; then #some files are already unzipped for some reason
		
		
		if [[ ! -e seq/${seq} ]]; then
			mkdir seq/$seq
		fi

		cp $aln_unzipped ${PWD}/seq/$seq
	fi
done


