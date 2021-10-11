#!/bin/bash
#
# Ingest SAdLSA data into the given database.
#
# SAdLSA output is generally stored in a set of directory hierarchies like this:
#
# {PROTEIN}/{SADLSA_RUN}/data_files
#
# E.g., WP_164928139.1/sadlsa_pdb70_210310/ will contain the SAdLSA score
# and alignment files for the sadlsa_pdb70_210310 SAdLSA run for protein
# WP_164928139.1.  The corresponding data files in that directory will be
# WP_164928139.1_aln.dat.gz and WP_164928139.1_sco.dat.gz.
#
# NOTE THAT THIS WILL APPEND DATA, SO IF THIS WAS ALREADY RUN, THEN THIS
# WILL JUST REDUNDANTLY ADD THE SAME DATA. Drop the tables before hand if you
# do not want to redundantly add this SAdLSA data.
# (TODO: explore use of UPSERTs to update data if it already exists, otherwise
# add it.)
#
# So we need to specify the directory of all the SAdLSA proteins and the
# name of the SAdLSA run.  Then we can grind through all the data to
# add it to the database.  We also need to know where the sqlite3 database
# is located to add to.
#
# -p protein directory
# -d sqlite3 database
# -s SAdLSA run directory

usage() {
  echo "ingest_sadlsa.sh -p <protein data directory> -d <sqlite3 database> -s <SAdLSA run directory>"
  exit 2
}

while getopts "hp:d:s:" arg; do
    case $arg in
    h)
      usage ;;
    p)
      PROTEIN_DIR=$OPTARG ;;
    d)
      DATABASE=$OPTARG ;;
    s)
      SADLSA_DIR=$OPTARG ;;
    esac
done

echo "PROTEIN_DIR is ${PROTEIN_DIR}"
echo "DATABASE is ${DATABASE}"
echo "SADLSA_DIR is ${SADLSA_DIR}"

# TODO use GNU parallel to signicantly speed this up
