#!/usr/bin/env python3
""" for ingesting TMalign data into sqlite3 database

TODO add support for alignment files; currently only ingests score files

"""
import argparse
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

import sys
from pathlib import Path
from time import time

from database.database import write_df_to_db
from parsers.tmalign_parser import parse_tmalign_score_file




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



    logging.info(f'Done ingesting TMalign data into {args.database} in '
                 f'{time() - start_time:0.2f} seconds')
