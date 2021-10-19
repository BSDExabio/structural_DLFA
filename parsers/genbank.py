#!/usr/bin/env python3
"""
    GenBank file support.
"""
from Bio import SeqIO


def protein_locus_dicts(genbank_file):
    """ This returns dicts for mapping between protein ID and locus tags from
        a Genbank file

    Note that we ignore the mapping from protein to locus tag for the old locus
    tags because then there would be multiple possible locus tags per protein.
    And, generally, we lookup locus tags for proteins, not the other way round.

    :param genbank_file: is the Genbank file from which we want to read
    :returns: two dicts, protein ID -> locus tag and locus tag -> protein ID
    """
    records = list(SeqIO.parse(genbank_file, "genbank"))

    locus_to_protein = {}
    protein_to_locus = {}

    for record in records:
        for feature in record.features:

            if 'protein_id' in feature.qualifiers:
                # We have a protein ID, but there may be multiple incarnations
                # for a corresponding locus tag in this record.  We try
                # all the possible combos and create corresponding look-up
                # entries for every single one of them.  Oh, and we strip out
                # any locus tag underscores so that we have a consistent,
                # universal standard for matching locus tags with their proteins
                if 'locus_tag' in feature.qualifiers:
                    for locus_tag in feature.qualifiers['locus_tag']:
                        locus_to_protein[locus_tag.replace('_','')] = \
                            feature.qualifiers['protein_id'][0]
                        protein_to_locus[feature.qualifiers['protein_id'][0]] = \
                            locus_tag.replace('_','')

                if 'old_locus_tag' in feature.qualifiers:
                    for old_locus_tag in feature.qualifiers['old_locus_tag']:
                        locus_to_protein[old_locus_tag.replace('_','')] = \
                            feature.qualifiers['protein_id'][0]
                        # This will *over-write* any previous mappings from
                        # the protein ID to locus tag, which we do not want.
                        # protein_to_locus[feature.qualifiers['protein_id'][0]] = \
                        #     old_locus_tag.replace('_','')
            elif 'inference' in feature.qualifiers:
                if 'locus_tag' in feature.qualifiers:
                    for locus_tag in feature.qualifiers['locus_tag']:
                        locus_to_protein[locus_tag.replace('_','')] = \
                            feature.qualifiers['inference'][0]
                        # Apparently this just means this is *similar* to a
                        # protein, but is not an actual mapping to a protein,
                        # therefore there is not a valid protein to locust tag
                        # mapping here.
                        # protein_to_locus[feature.qualifiers['sequence:RefSeq:'][0]] = \
                        #     locus_tag.replace('_','')

                if 'old_locus_tag' in feature.qualifiers:
                    for old_locus_tag in feature.qualifiers['old_locus_tag']:
                        locus_to_protein[old_locus_tag.replace('_','')] = \
                            feature.qualifiers['inference'][0]


    return protein_to_locus, locus_to_protein


if __name__ == '__main__':
    # Test harness for Genebank parsing

    protein_to_locus, locus_to_protein = protein_locus_dicts(
        '/Users/may/Downloads/alphafold/GCF_000195755.1_ASM19575v1_genomic.gbff')

    assert protein_to_locus['WP_010937312.1'] == 'DVURS00005'
    assert locus_to_protein['DVURS00005'] == 'WP_010937312.1'

    # Catch stupid edge case
    assert locus_to_protein['DVURS05145'] == 'COORDINATES: similar to AA sequence:RefSeq:WP_011792536.1'
    assert locus_to_protein['DVU1087'] == 'COORDINATES: similar to AA sequence:RefSeq:WP_011792536.1'
