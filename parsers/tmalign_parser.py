#!/usr/bin/env python3
"""
    For parsing TMalign-Apoc data.  It expects data like this:

tname  malnlen  mrmsd   mscore  mseqid	Description
3LFZ_A     125  1.700 0.76562  28.8%	| protein MJ1225; MJ1225, AMPK, AMP, ADP, CBS; HET: ATP, ADP, AMP; 2.2A {Methanocaldococcus jannaschii}
3KH5_A     125  1.710 0.76354  28.8%	| protein MJ1225; MJ1225, AMPK, AMP, ADP, CBS; HET: ADP, AMP; 2.1A {Methanocaldococcus jannaschii}
"""
from pathlib import Path
import pandas as pd
import gzip

def _get_protein(fn):
    """ return the protein name, which is part of the filename
        E.g. 'WP_164928147.1_aln.dat' -> 'WP_164928147.1'

        TOTO consolidate this into a utils.py since this function is used in
        more than one parser.
    """
    # Just work backwards to the first '_' in just the file name, and there's
    # your snipping point
    name = Path(fn).name  # We do this to snip any leading directory paths
    return name[:name.rfind('_')]


def parse_tmalign_score_file(fn):
    """ Read given TMalign scores file into a pandas dataframe

    :param fn: of score file
    :return: pandas dataframe
    """
    columns = ['protein', 'tname', 'malnlen', 'mrmsd', 'mscore', 'mseqid', 'Description']

    rows = [] # will contain rows of dicts corresponding to score data

    protein = _get_protein(fn)

    if fn[-1] == 'z':
        # compressed file names should end in 'z'
        with gzip.open(fn, mode='rt') as f:
            lines = f.readlines()
    else:
        # Someone manually uncompressed the file, probably to look at it
        with open(fn, mode='rt') as f:
            lines = f.readlines()

    for line in lines[1:]: # 1: to skip header

        temp = line.strip().split('|')
        rec = [protein] + temp[0].strip().split() + [temp[1].strip()]   # prepping the list of alignment result elements

        # This is just a fast way of mapping the column names to each
        # successive list element.  We could have laboriously assigned each
        # dictionary element, instead, but why do that?  :P
        temp_dict = {index: value for index, value in zip(columns, rec)}

        # Drop the noisome '%' so that later we can do numeric operations on it
        temp_dict['mseqid'] = temp_dict['mseqid'][:-1]
        rows.append(temp_dict.copy())

    df = pd.DataFrame(rows)

    # Converts from objects to strings
    df = df.convert_dtypes()

    # Convert to numeric those things that are numeric
    df[[ 'malnlen' , 'mrmsd'  , 'mscore',  'mseqid']] = \
        df[['malnlen', 'mrmsd', 'mscore', 'mseqid']].apply(pd.to_numeric)

    return df



def parse_tmalign_align_file(fn):
    """ Read given TMalign alignment file into a pandas dataframe

    :param fn: of score file
    :return: pandas dataframe
    """
    # TODO implement
    pass


if __name__ == '__main__':
    # test harness for parsers

    df = parse_tmalign_score_file('/Users/may/Projects/data/Alphafold/TMalign-APoc/WP_010939197.1/WP_010939197.1_sco.dat')

    pass
