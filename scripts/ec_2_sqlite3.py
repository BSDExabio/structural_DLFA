#!/usr/bin/env python3
"""
    Imports ENZYME nomenclature database to an sqlite3 table

usage: ec_2_sqlite3.py [-h] --ec-file EC_FILE --database DATABASE

ENZYME nomenclature database importer

optional arguments:
  -h, --help           show this help message and exit
  --ec-file EC_FILE    ENZYME nomenclature database file
  --database DATABASE  sqlite3 in which to insert EC data
"""
import argparse
import logging
import string
import sys
from pathlib import Path

import pandas as pd
from rich import pretty
pretty.install()
from rich.logging import RichHandler

from rich.traceback import install
install()

from database.database import write_df_to_db

rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='INFO', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])
logger = logging.getLogger(__name__)

# name of target database table
EC_TABLE = 'enzyme_desciptions'


def ec_file_to_df(ec_file):
    """ parse the file containing the ENZYME nomenclature database

    Outside of the header and footer of this plain text file, the data looks
    like this:

    1. -. -.-  Oxidoreductases.
    1. 1. -.-   Acting on the CH-OH group of donors.
    1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.

    So: columns 0 is class, 2-3 subclass, 5-6 sub-subclass, 11+: is description

    :param ec_file: EC file of nomenclature database
    :return: data frame of EC file
    """
    rows = []  # will contain rows of dicts corresponding to EC data_frame

    # Gets updated for each entry, then appended to rows
    curr_entry = {'class'           : 0,
                  'subclass'        : 0,
                  'subsubclass'     : 0,
                  'class_name'      : '',
                  'subclass_desc'   : '',
                  'subsubclass_desc': ''}

    with ec_file.open('r') as ec_db:
        for line in ec_db:
            if line[0] in string.digits:
                # we have a hit on a line of data
                if line[6] in string.digits:
                    # this is a sub-sub-class
                    curr_entry['subsubclass'] = line[5:7]
                    curr_entry['subsubclass_desc'] = line[11:]
                    rows.append(curr_entry.copy())
                elif line[3] in string.digits:
                    # this is a sub-class
                    curr_entry['subclass'] = line[2:4]
                    curr_entry['subclass_desc'] = line[11:]
                else:
                    # this is a class
                    curr_entry['class'] = line[0]
                    curr_entry['class_name'] = line[11:]

    ec_df = pd.DataFrame(rows)

    # Strip out whitespace
    # (https://stackoverflow.com/questions/33788913/pythonic-efficient-way-to-strip-whitespace-from-every-pandas-data-frame-cell-tha/33789292)
    ec_df = ec_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # Converts from objects to strings
    ec_df = ec_df.convert_dtypes()

    # But these are really integers, and not very big ones, so Int8
    # is big enough
    ec_df['class'] = pd.to_numeric(ec_df['class']).astype('Int8')
    ec_df['subclass'] = pd.to_numeric(ec_df['subclass']).astype('Int8')
    ec_df['subsubclass'] = pd.to_numeric(ec_df['subsubclass']).astype('Int8')

    # Snip out those awful trailing periods
    ec_df['class_name'] = ec_df['class_name'].str.replace(r'.$','')
    ec_df['subclass_desc'] = ec_df['subclass_desc'].str.replace(r'.$','')
    ec_df['subsubclass_desc'] = ec_df['subsubclass_desc'].str.replace(r'.$','')

    return ec_df





if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='ENZYME nomenclature database importer')
    parser.add_argument('--ec-file', required=True,
                        help='ENZYME nomenclature database file')
    parser.add_argument('--database', required=True,
                        help='sqlite3 in which to insert EC data')

    args = parser.parse_args()

    ec_file = Path(args.ec_file)
    if not ec_file.exists():
        logging.critical(f'{args.ec_file} does not exist ... exiting')
        sys.exit(1)

    logging.info(f'Ingesting {args.ec_file}')

    ec_df = ec_file_to_df(ec_file)

    write_df_to_db(ec_df, EC_TABLE, args.database)

    logging.info(f'Done.  Added {len(ec_df)} entries to'
                 f' {args.database} table {EC_TABLE}.')
