#!/usr/bin/env python3
"""
    For parsing TMalign data.

"""
from pathlib import Path
import pandas as pd
import gzip

def _get_protein(fn):
    """ return the protein name, which is part of the filename path
        Unfortunately, there are now multiple formats for filename paths.
        The function needs to find the protein name for all formats.
            REALLY OLD: protein name in fn's file name
                '~/WP_164928147.1_aln.dat' -> 'WP_164928147.1'
            OLD: protein name in fn's parent directory name
                '~/WP_014320981.1/ranked_alignment_results.dat' -> 'WP_014320981.1'
            NEW: protein name in fn's grandparent directory name
                '~/WP_014320981.1/pdb70_2022_03_19_rbd/ranked_alignment_results.dat' -> 'WP_014320981.1'
    """
    return Path(fn).parent.parent.name
    #return Path(fn).parent.name


def parse_tmalign_score_file(fn):
    """ Read given TMalign scores file into a pandas dataframe
    These files are written in csv format to begin with so will
    be easy to read into a pandas dataframe

    :param fn: path to the score file
    :return: pandas dataframe
    """
    
    protein = _get_protein(fn)

    df = pd.read_csv(fn)
    df['protein']    = [protein for i in range(len(df))]
    df['Description'] = ['' for i in range(len(df))]
    
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

    df = parse_tmalign_score_file('/Users/may/Projects/data/Alphafold/TMalign/WP_010939197.1/TMalign_cp_pdb70_2022_03_19/ranked_alignment_results.dat')

    pass
