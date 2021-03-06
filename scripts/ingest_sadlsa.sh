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
# -c number of cores to use in parallel

usage() {
  echo "ingest_sadlsa.sh -p <protein data directory> -d <sqlite3 database> -s <SAdLSA run directory>"
  exit 2
}

# Default to 8 cores; parallel normally figures out the correct number of
# cores, but not on Macs, apparently.  :P
CORES=8

while getopts "hp:d:s:" arg; do
    case $arg in
    h)
      usage ;;
    p)
      PROTEIN_DIR=$OPTARG ;;
    d)
      DATABASE=$OPTARG ;;
    c)
      CORES=$OPTARG ;;
    s)
      SADLSA_DIR=$OPTARG ;;
    esac
done

echo "PROTEIN_DIR is ${PROTEIN_DIR}"
echo "DATABASE is ${DATABASE}"
echo "SADLSA_DIR is ${SADLSA_DIR}"
echo "CORES is ${CORES}"

# So we can find database import
export PYTHONPATH=..:../database:../parsers:$PYTHONPATH

# Replaced by GNU parallel to significantly speed this up
#for f in ${PROTEIN_DIR}/*/${SADLSA_DIR}/ ; do
#  python3 ./sadlsa_2_sqlite3.py --database $DATABASE --sadlsa-dir $f
#done

# We need to run the first job solo so that any sqlite3 database tables are
# created.  Then we can run parallel on the rest.
ls -d ${PROTEIN_DIR}/*/${SADLSA_DIR}/ > /tmp/FILES

# Just do the first to build the tables; if we didn't do this and tried to do
# everything in a single parallel call, the first $CORE processes would try
# to build any new tables a the same time, which is bad.  So just let one
# create any needed database tables so that we can then later run parallel
# safely since then all those tables are guaranteed to exist.
python3 ./sadlsa_2_sqlite3.py --database $DATABASE --sadlsa-dir `head -1 /tmp/FILES`

# Then crank on all the rest with no worries of creating tables, which were done
# in the previous line.
tail -n +2 /tmp/FILES | parallel -j $CORES python3 ./sadlsa_2_sqlite3.py --database $DATABASE --sadlsa-dir {1}
