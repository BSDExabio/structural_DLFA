#!/usr/bin/env python3
"""
    Imports ENZYME nomenclature database to an sqlite3 table
"""
import argparse
import re
import string
from pathlib import Path
import pandas as pd
import numpy as np
import sys
from sqlalchemy import create_engine
import logging
from rich.logging import RichHandler
from rich.table import Table
from rich import print
from rich import pretty

pretty.install()

from rich.traceback import install

install()

rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='DEBUG', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])
logger = logging.getLogger(__name__)


def ec_file_to_df(ec_file):
    """ parse the file containing the ENZYME nomenclature database

    Outside of the header and footer of this plain text file, the data looks
    like this:

    1. -. -.-  Oxidoreductases.
    1. 1. -.-   Acting on the CH-OH group of donors.
    1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.

    So: columns 0 is class, 2-3 subclass, 5-6 sub-subclass, 11: is description

    :param ec_file: EC file of nomenclature database
    :return: data frame of EC file
    """
    ec_df = pd.DataFrame(columns=['class', 'subclass', 'subsubclass',
                                  'class_name', 'sub_class_desc',
                                  'sub_sub_class_desc'])

    curr_class_idx = 0 # index for current class
    curr_subclass_idx = 0 # index for current subclass
    curr_class_name = '' # name of current class
    curr_subclass_desc = '' # current subclass description

    with ec_file.open('r') as ec_db:
        while line in ec_db:
            if line[0] in string.digits:
                # we have a hit on a line of data
                if line[6] in string.digits:
                    # this is a sub-sub-class
                    pass
                elif line[3] in string.digits:
                    # this is a sub-class
                    pass
                else:
                    # this is a class
                    curr_class_idx = int(line[0])
                    ec_df[curr_subclass] = {'class' : line[11:]}



    return ec_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ENZYME nomenclature database importer')
    parser.add_argument('--ec-file',  required=True,
                        help='ENZYME nomenclature database file')
    parser.add_argument('--database', required=True,
                        help='sqlite3 in which to insert EC data')

    args = parser.parse_args()

    ec_file = Path(args.ec_file)

    if not ec_file.exists():
        print(f'{args.ec_file} does not exist ... exiting')
        sys.exit(1)

    entries = ec_file_to_df(ec_file)
