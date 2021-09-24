#!/usr/bin/env python3
"""
    Imports ENZYME nomenclature database to an sqlite3 table
"""
import argparse
import re
from pathlib import Path
import pandas as pd
import sys
from sqlalchemy import create_engine
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


def parse_ec_file(ec_file):
    """ parse the file containing the ENZYME nomenclature database

    :param ec_file: EC file of nomenclature database
    :return: dict of nomenclature hierarchy
    """
    entries = {}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='ENZYME nomenclature database importer')
    parser.add_argument('--ec-file', help='ENZYME nomenclature database file')
    parser.add_argument('--database', help='sqlite3 in which to insert EC data')

    args = parser.parse_args()

    ec_file = Path(args.ec_file)

    if not ec_file.exists():
        print(f'{args.ec_file} does not exist ... exiting')
        sys.exit(1)

    entries = parse_ec_file(ec_file)
