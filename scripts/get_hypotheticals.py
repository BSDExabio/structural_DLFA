#!/usr/bin/env python3
""" This will extract any hypothetical proteins from a given Genbank file.

"""
import argparse
import logging
import string
import sys
import re
from pathlib import Path

from rich import pretty

pretty.install()

from rich.logging import RichHandler

rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='INFO', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])
logger = logging.getLogger(__name__)

from rich.traceback import install

install()

from Bio import SeqIO

DESCRIPTION = """Find hypothetical proteins in Genbank files

Find hypothetical proteins in Genbank files and either write them to a 
corresponding FASTA file, or echo the protein IDs to stdout.  I.e., 
if a FASTA output file name is not given, then just print the protein 
IDs to stdout. 
"""


def find_hypothetical_proteins(genbank_file):
    """ Return any hypothetical proteins that are NOT psuedo-genes in the given
        Genbank file.

    :param genbank_file: in which we'll be looking for hypothetical proteins
    :return: list of matching records
    """
    records = list(SeqIO.parse(genbank_file, "genbank"))

    matches = []

    for record in records:
        for feature in record.features:

            if 'protein_id' in feature.qualifiers and \
                    'product' in feature.qualifiers:
                # We look for the work "hypothetical" anywhere in the product
                # description.  Note this implicitly filters out pseudo-genes
                # since those do not have a `protein_id`.
                if re.search(r'hypothetical', feature.qualifiers['product'][0]):
                    matches.append(feature)
    return matches


def print_matches(matches):
    """ Blat out matches to stdout """
    for match in matches:
        print(match.qualifiers['protein_id'][0])


def write_to_fasta(matches, fasta_filename):
    """ Write matches to FASTA file """
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('genbank_file',
                        help='Genbank file in which we want to find '
                             'hypothetical proteins')
    parser.add_argument('--fasta-file',
                        help='Where we optionally want to write '
                             'hypothetical proteins in FASTA format')

    args = parser.parse_args()

    genbank_file = Path(args.genbank_file)

    if not genbank_file.exists():
        logger.critical(f'{args.genbank_file} does not exist')

    matches = find_hypothetical_proteins(genbank_file)

    if matches == []:
        logger.warning(f'No hypothetical proteins found in {args.genbank_file}')
        sys.exit(1)

    if args.fasta_file is None:
        # Just blat out the protein IDs to stdout since we didn't specify
        # a FASTA file
        print_matches(matches)
    else:
        # TODO implement writing to FASTA file
        pass
