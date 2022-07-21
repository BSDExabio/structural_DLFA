#!/bin/bash
#
# Ingest TMalign data into the given database.
#
# TMalign output is generally stored in a set of directory hierarchies like this:
#
# run/{PROTEIN}/data_files
#
# E.g., WP_010939197.1 will contain the TMalign score
# and alignment files for the TMalign_pdb70_210310 TMalign run for protein
# WP_010939197.1.  The corresponding data files in that directory will be
# WP_010939197.1_aln.dat and WP_010939197.1_sco.dat.
#
# NOTE THAT THIS WILL APPEND DATA, SO IF THIS WAS ALREADY RUN, THEN THIS
# WILL JUST REDUNDANTLY ADD THE SAME DATA. Drop the tables before hand if you
# do not want to redundantly add this TMalign data.
# (TODO: explore use of UPSERTs to update data if it already exists, otherwise
# add it.)
#
# So we need to specify the directory of all the TMalign proteins and the
# name of the TMalign run.  Then we can grind through all the data to
# add it to the database.  We also need to know where the sqlite3 database
# is located to add to.
#
# -d sqlite3 database
# -t TMalign run directory
# -c number of cores to use in parallel

usage() {
  echo "ingest_tmalign.sh -p <protein data directory> -d <sqlite3 database> -s <TMalign run directory>"
  exit 2
}

# Default to 8 cores; parallel normally figures out the correct number of
# cores, but not on Macs, apparently.  :P
CORES=8

while getopts "hp:d:s:" arg; do
    case $arg in
    h)
      usage ;;
    d)
      DATABASE=$OPTARG ;;
    c)
      CORES=$OPTARG ;;
    t)
      TMalign_DIR=$OPTARG ;;
    esac
done

echo "DATABASE is ${DATABASE}"
echo "TMalign_DIR is ${TMalign_DIR}"
echo "CORES is ${CORES}"

# So we can find database import
export PYTHONPATH=..:../database:../parsers:$PYTHONPATH

# We need to run the first job solo so that any sqlite3 database tables are
# created.  Then we can run parallel on the rest.
ls -d ${TMalign_DIR}/WP* > /tmp/FILES

# Just do the first to build the tables; if we didn't do this and tried to do
# everything in a single parallel call, the first $CORE processes would try
# to build any new tables a the same time, which is bad.  So just let one
# create any needed database tables so that we can then later run parallel
# safely since then all those tables are guaranteed to exist.
python3 ./tmlign_2_sqlite3.py --database $DATABASE --tmalign-dir `head -1 /tmp/FILES`

# Then crank on all the rest with no worries of creating tables, which were done
# in the previous line.
tail -n +2 /tmp/FILES | parallel -j $CORES python3 ./tmalign_2_sqlite3.py --database $DATABASE --tmalign-dir {1}
