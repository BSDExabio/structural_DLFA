#!/usr/bin/env python3
""" for ingesting TMalign data into sqlite3 database

TODO add support for alignment files; currently only ingests score files

usage: tmalign_2_sqlite3.py [-h] --tmalign-dir TMALIGN_DIR --database DATABASE

TMalign data importer

optional arguments:
  -h, --help            show this help message and exit
  --tmalign-dir TMALIGN_DIR
                        directory containing TMalign alignments and scores
  --database DATABASE   sqlite3 in which to insert TMalign data
"""
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

import sys
from pathlib import Path
from time import time

from database.database import write_df_to_db
from parsers.tmalign_parser import parse_tmalign_score_file


# sqlite3 table that will contain all the scores
TMALIGN_SCORES_TABLE = 'tmalign_scores'


def add_scores_to_db(tmalign_dir, database):
    """ Add TMalign score to database for given TMalign data directory

    The data directory should contain two files, one for scores, and the other
    for alignments.  We want the scores file, which should match the pattern
    "*_sco.dat8*".

    :param tmalign_dir: output directory containing TMalign scores files
    :param database: path to sqlite3 database to write the scores table
    :return: None
    """
    # grab score file first
    score_file = list(sadlsa_dir.glob('*_sco.dat*'))
    if [] == score_file:
        logging.critical(
            f'No score file found in path {sadlsa_dir!s} ... exiting')
        sys.exit(1)

    # We return the first element of the glob because there should only be one
    # score file.
    logging.info(f'Reading {score_file[0]}')
    tmalign_score_df = parse_tmalign_score_file(str(score_file[0]))

    write_df_to_db(data_frame=tmalign_score_df,
                   table=TMALIGN_SCORES_TABLE,
                   database=database)

    logging.info(f'Added {len(tmalign_score_df)} entries to'
                 f' {database} table {TMALIGN_SCORES_TABLE}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='TMalign data importer')
    parser.add_argument('--tmalign-dir', required=True,
                        help='directory containing TMalign alignments and scores')
    parser.add_argument('--database', required=True,
                        help='sqlite3 in which to insert TMalign data')

    args = parser.parse_args()

    tmalign_dir = Path(args.tmalign_dir)
    if not tmalign_dir.exists():
        logging.critical(f'{args.tmalign_dir} does not exist ... exiting')
        sys.exit(1)

    logging.info(f'Ingesting {args.tmalign_dir}')
    start_time = time()

    add_scores_to_db(tmalign_dir, args.database)

    # TODO add alignments to database

    logging.info(f'Done ingesting TMalign data into {args.database} in '
                 f'{time() - start_time:0.2f} seconds')
