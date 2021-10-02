#!/usr/bin/env python3
"""
    For importing SAdLSA data into pandas dataframes.

    This is used by another utility for importing into sqlite3 database.
"""
from pathlib import Path
import pandas as pd


def _get_protein(fn):
    """ return the protein name, which is part of the filename

        E.g. 'WP_164928147.1_aln.dat' -> 'WP_164928147.1'
    """
    # Just work backwards to the first '_' in just the file name, and there's
    # your snipping point
    name = Path(fn).name # We do this to snip any leading directory paths
    return name[:name.rfind('_')]


def parse_sadlsa_score_file(fn, count=10):
    """ translate given SAdLSA score file into pandas dataframe

    :param fn: string filename of score file
    :param count: how many of the top scores we want to get
    :return: pandas dataframe of
    """
    import gzip

    columns = ["protein", "seqname", "alnlen", "range1", "tmscore1", "range2",
               "tmscore2", "alnscore", "seq_id", "desc"]
    df = pd.DataFrame(columns=columns)

    # extract the protein name from the filename
    protein = _get_protein(fn)

    if fn[-1] == 'z':
        # compressed file names should end in 'z'
        with gzip.open(fn, mode='rt') as f:
            lines = f.readlines()
    else:
        # Someone manually uncompressed the file, probably to look at it
        with open(fn, mode='rt') as f:
            lines = f.readlines()

    for line in lines[3:count + 3]:
        temp = line.strip().split('\t|')
        rec = [protein] + temp[0].split()
        rec.append(temp[1])

        series = pd.Series(rec, index=df.columns)
        df = df.append(series, ignore_index=True)

    return df


def parse_sadlsa_aln_file(fn, count=10):
    """ translate given SAdLSA alignment file into pandas dataframe

    TODO consider splitting this into two dataframes, one for the alignment
    block header, and the other for the rows of alignment data; that should
    save redundant column information.

    :param fn: string filename of alignment file
    :param count: how many of the top scores we want to get
    :return: pandas dataframe of
    """
    import gzip
    import re

    prog = re.compile(
        r'### Alignment (\d+) to: (......) naln=(\d+) score=(\d+\.\d*) tms1=(\d+\.\d*) tms2=(\d+\.\d*) sid=(\d+\.\d*)')

    columns = ["protein", "Aln_num", "Prot_ID", "naln", "score", "tms1", "tms2",
               "sid", "Ind", "Res1", "AA1", "Res2", "AA2", "MeanDist", "Bin",
               "Prob<3", "Prob<5", "Prob<8"]
    df = pd.DataFrame(columns=columns)

    # extract the protein name from the filename
    protein = _get_protein(fn)

    if fn[-1] == 'z':
        # compressed file names should end in 'z'
        with gzip.open(fn) as f:
            lines = f.readlines()
    else:
        # Someone manually uncompressed the file, probably to look at it
        with open(fn, mode='rt') as f:
            lines = f.readlines()

    for line in lines:

        search_results = prog.search(line)
        if search_results:
            # We have a hit on the start of an alignment block
            if len(df) == count:
                # We're done if we have the user's target number of rows
                break
            alignment_header_data = list(search_results.groups())
        elif '#Ind' == line[0:4]:
            # This is the header for the alignment block, so we can skip it
            continue
        elif line.strip() == '':
            # This is a line between blocks, so we can just skip it
            continue
        else:
            # We have a row of actual data, so split out the data (after
            # stripping out any '*' in the line) and append that to the
            # alignment block info
            rec = [protein] + alignment_header_data + \
                  line.strip().replace('*', '').split()
            series = pd.Series(rec, index=df.columns)
            df = df.append(series, ignore_index=True)

    return df


if __name__ == '__main__':
    # test harness for these functions

    # Where the score and alignment files of interest are.
    base_path = Path('/Users/may/Projects/data/PSP/desulfovibrio_vulgaris/out/'
                     'WP_164928147.1/sadlsa_pdb70_210310')

    print('Reading score file')
    score_df = parse_sadlsa_score_file(
        str(base_path / 'WP_164928147.1_sco.dat'))

    print('Reading alignment file')
    align_df = parse_sadlsa_aln_file(str(base_path / 'WP_164928147.1_aln.dat'))

    pass
