#!/usr/bin/env
"""
    This script will take a file that contains one or more columns of locus tags
    and will add columns for the corresponding protein IDs.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

from parsers.genbank import protein_locus_dicts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Look up protein IDs from locus tags')
    parser.add_argument('--infile', '-f', required=True,
                        help='File contains locus tags to be looked up')
    parser.add_argument('--genbank', '-g', required=True,
                        help='Genbank file that contain protein IDs and locust tags')
    parser.add_argument('tag-columns', nargs='+',
                        help='Identify what columns are for locust tags')

    args = parser.parse_args()

    infile = Path(args.infile)
    genbank_file = Path(args.genbank)

    if not infile.exists():
        print(f'{args.infile} does not exist ... exiting')
        sys.exit(1)

    if not genbank_file.exists():
        print(f'{args.genbank} does not exist ... exiting')
        sys.exit(1)

    data = pd.read_csv(str(infile), delim_whitespace=True)

    protein_to_locus, locus_to_protein = protein_locus_dicts(str(genbank_file))

    pass
