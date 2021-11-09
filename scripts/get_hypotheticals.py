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


def _is_hypothetical_protein(feature):
    """ returns True if hypothetical protein """
    return 'protein_id' in feature.qualifiers and \
           'product' in feature.qualifiers and \
           re.search(r'hypothetical', feature.qualifiers['product'][0])


def write_stdout(genbank_file):
    """ Blat out matches to stdout """

    def print_id(feature):
        print(feature.qualifiers['protein_id'][0])

    records = SeqIO.parse(genbank_file, "genbank")

    for record in records:
        [print_id(feature) for feature in record.features
         if _is_hypothetical_protein(feature)]


def write_to_fasta(genbank_file, fasta_filename):
    """ Write matches to FASTA file """
    # First pull everything into a dictionary, and then splice out all the
    # feature records that DON'T have hypothetical proteins.  We'll later
    # write that dictionary to a FASTA file.
    records = SeqIO.to_dict(SeqIO.parse(genbank_file, "genbank"))

    for record in records:
        # Filter for features that just have hypothetical proteins that are
        # NOT pseudo-genes.
        new_features = [feature for feature in record.features
                        if _is_hypothetical_protein(feature)]
        record.features = new_features

    with open(fasta_filename, 'w') as fasta_file:
        count = SeqIO.write(records, fasta_file, 'fasta')

    logger.info(f'Wrote {count} records to {fasta_filename}')


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

    # matches = find_hypothetical_proteins(genbank_file)

    # if matches == []:
    #     logger.warning(f'No hypothetical proteins found in {args.genbank_file}')
    #     sys.exit(1)

    if args.fasta_file is None:
        # Just blat out the protein IDs to stdout since we didn't specify
        # a FASTA file
        write_stdout(genbank_file)
    else:
        write_to_fasta(genbank_file, args.fasta_file)
