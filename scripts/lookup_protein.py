#!/usr/bin/env
"""
    This script will take a file that contains one or more columns of locus tags
    and will add columns for the corresponding protein IDs.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Look up protein IDs from locus tags')
    parser.add_argument('infile', help='File contains locus tags to be looked up')

    args = parser.parse_args()

    infile = Path(args.infile)

    if not infile.exists():
        print(f'{args.infile} does not exist ... exiting')
        sys.exit(1)

    data = pd.read_csv(str(infile), delim_whitespace=True)

    pass
