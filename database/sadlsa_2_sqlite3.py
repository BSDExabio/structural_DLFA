#!/usr/bin/env python3
"""
    Imports SAdLSA data into an sqlite3 table
"""
import argparse
import logging
import sys
from pathlib import Path

from database import write_df_to_db
from parsers.sadlsa_parsers import parse_sadlsa_score_file, parse_sadlsa_aln_file

# sqlite3 table that will contain all the scores
SADLSA_SCORES_TABLE = 'sadlsa_scores'

# sqlite3 table that will contain all the alignment data
SADLSA_ALIGNMENTS_TABLE = 'sadlsa_alignments'


def add_scores_to_db(sadlsa_dir, database):
    """ Add SAdLSA score to database for given SAdLSA output directory

    The output directory should contain two files, one for scores, and the other
    for alignments.  We want the scores file, which should match the pattern
    "*_sco.dat8*".

    :param sadlsa_dir: output directory containing scores file
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
    sadlsa_score_df = parse_sadlsa_score_file(str(score_file[0]))
    write_df_to_db(data_frame=sadlsa_score_df,
                   table=SADLSA_SCORES_TABLE,
                   database=database)

    logging.info(f'Done.  Added {len(sadlsa_score_df)} entries to'
                 f' {database} table {SADLSA_SCORES_TABLE}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='SAdLSA data importer')
    parser.add_argument('--sadlsa-dir', required=True,
                        help='directory containing SAdLSA alignment and score filed got a single protein')
    parser.add_argument('--database', required=True,
                        help='sqlite3 in which to insert EC data')

    args = parser.parse_args()

    sadlsa_dir = Path(args.sadlsa_dir)
    if not sadlsa_dir.exists():
        logging.critical(f'{args.sadlsa_dir} does not exist ... exiting')
        sys.exit(1)

    logging.info(f'Ingesting {args.sadlsa_dir}')

    add_scores_to_db(sadlsa_dir, args.database)
