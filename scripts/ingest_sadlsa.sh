#!/usr/bin/env bash
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
# Call getopt to validate the provided input.
options=$(getopt -o hpds: -- "$@")
[ $? -eq 0 ] || {
    echo "Incorrect options provided"
    exit 1
}
eval set -- "$options"
while true; do
    case "$1" in
    -h)
      echo "Usage: ingest_sadlsa.sh -p <PROTEIN> -d <DATABASE> -s <SADLSA_DIR>"
      exit 1
      ;;
    -p)
      PROTEIN_DIR=$1;
      ;;
    -d)
      DATABASE=$1;
      ;;
    -s)
      SADLSA_DIR=$1;
      ;;
    --)
      shift;
      break
      ;;
    esac
    shift
done

if [ ${PROTEIN_DIR}x = x ]; then
  echo "Missing protein dir option -p"
  exit
fi
if [ ${DATABASE}x = x ]; then
  echo "Missing protein dir option -d"
  exit
fi
if [ ${SADLSA_DIR}x = x ]; then
  echo "Missing protein dir option -s"
  exit
fi

echo "PROTEIN_DIR is ${PROTEIN_DIR}"
echo "DATABASE is ${DATABSE}"
echo "SADLSA_DIR is ${SADLSA_DIR}"
