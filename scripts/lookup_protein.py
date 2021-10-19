#!/usr/bin/env python3
"""
    This script will take a file that contains one or more columns of locus tags
    and will add columns for the corresponding protein IDs.

    usage: lookup_protein.py [-h] --infile INFILE --genbank GENBANK tag_columns [tag_columns ...]

    Look up protein IDs from locus tags

    positional arguments:
      tag_columns           Identify what columns are for locust tags

    optional arguments:
      -h, --help            show this help message and exit
      --infile INFILE, -f INFILE
                            File contains locus tags to be looked up
      --genbank GENBANK, -g GENBANK
                            Genbank file that contain protein IDs and locust tags
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
    parser.add_argument('tag_columns', nargs='+',
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

    df = pd.read_csv(str(infile), delim_whitespace=True)

    # Get two dicts, one for looking up locus tags, and the other for looking
    # up protein IDs.  We're only going to use the latter for this script,
    # though.
    protein_to_locus, locus_to_protein = protein_locus_dicts(str(genbank_file))

    # Now add columns for each locus tag column
    for locus_tag in args.tag_columns:
        protein_column = locus_tag + '_protein'

        df[protein_column] = \
            df[locus_tag].apply(lambda x: locus_to_protein.setdefault(x, pd.NA))

    print(df.to_csv(na_rep='No match', index=False))
