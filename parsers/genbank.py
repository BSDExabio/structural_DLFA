#!/usr/bin/env python3
"""
    GenBank file support.
"""
from Bio import SeqIO


def protein_locus_dicts(infile):
    """ This returns dicts for mapping between protein ID and locus tags from
        a Genbank file

    :param infile: is the Genbank file from which we want to read
    :returns: two dicts, protein ID -> locus tag and locus tag -> protein ID
    """
    records = SeqIO.parse(infile, "genbank")

    locus_to_protein = {}
    protein_to_locus = {}

    for record in records:
        for feature in record.features:

            if 'protein_id' in feature.qualifiers and \
                    'locus_tag' in feature.qualifiers:
                # Update the dicts with both the protein ID and locus tag; we
                # also strip out the bothersome underscore from the locus tag.
                locus_to_protein[feature.qualifiers['locus_tag'][0].replace('_','')] = \
                    feature.qualifiers['protein_id'][0]
                protein_to_locus[feature.qualifiers['protein_id'][0]] = \
                    feature.qualifiers['locus_tag'][0].replace('_','')

    return protein_to_locus, locus_to_protein


if __name__ == '__main__':
    # Test harness for Genebank parsing

    protein_to_locus, locus_to_protein = protein_locus_dicts(
        '/Users/may/Downloads/alphafold/GCF_000195755.1_ASM19575v1_genomic.gbff')

    assert protein_to_locus['WP_010937312.1'] == 'DVURS00005'
    assert locus_to_protein['DVURS00005'] == 'WP_010937312.1'
