#!/usr/bin/env python3
"""
    For importing APOC data into pandas dataframes.
    This is used by another utility for importing into sqlite3 database.
"""
import pandas as pd
from pathlib import Path

def _get_protein(fn):
    """ return the protein name, which is part of the filename
        E.g. 'WP_164928147.1_aln.dat' -> 'WP_164928147.1'
    """
    # Just work backwards to the first '_' in just the file name, and there's
    # your snipping point
    name = Path(fn).name  # We do this to snip any leading directory paths
    return name[:name.rfind('_')]


def parse_apoc_score_file(fn, count=10):
    """ translate given APOC score file into pandas dataframe
    :param fn: string filename of score file
    :param count: how many of the top scores we want to get
    :return: pandas dataframe of score results for top count results
    """
    import gzip

    columns = ["protein", "seqname", "tmalnlen", "tmrmsd", "tmscore", "seq_id", "desc"]

    rows = [] # will contain rows of dicts corresponding to score data

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

    for line in lines[1:count + 1]:
        temp = line.strip().split('\t|')
        rec = [protein] + temp[0].split() + [temp[1]]   # prepping the list of alignment result elements

        # This is just a fast way of mapping the column names to each
        # successive list element.  We could have laboriously assigned each
        # dictionary element, instead, but why do that?  :P
        temp_dict = {index: value for index, value in zip(columns, rec)}

        # Drop the noisome '%' so that later we can do numeric operations on it
        temp_dict['seq_id'] = temp_dict['seq_id'][:-1]
        rows.append(temp_dict.copy())

    df = pd.DataFrame(rows)

    # Converts from objects to strings
    df = df.convert_dtypes()

    # Convert to numeric those things that are numeric
    df[["tmalnlen", "tmrmsd", "tmscore", "seq_id"]] = df[["tmalnlen", "tmrmsd", "tmscore", "seq_id"]].apply(pd.to_numeric)

    return df


def parse_apoc_aln_file(fn, count=10):
    """ translate given APOC alignment file into pandas dataframe
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
        r'### Alignment (\d+) to: (\w+) naln=(\d+) score=(\d+\.\d*) sid=(\d+\.\d*)')

    # extract the protein name from the filename
    protein = _get_protein(fn)

    rows = []  # will contain rows of dicts corresponding to alignment data

    # Gets updated for each entry, then appended to rows
    curr_entry = {"protein" : protein,
                  "Aln_num" : 0,
                  "Prot_ID" : "",
                  "naln"    : 0,
                  "score"   : 0.0,
                  "sid"     : "",
                  "Ind"     : 0,
                  "Ch1"     : 0,
                  "Res1"    : 0,
                  "AA1"     : 0,
                  "Ch2"     : 0,
                  "Res2"    : 0,
                  "AA2"     : 0,
                  "Dist"    : 0,
                  "Cos"     : 0,
                  "S"       : 0}

    if fn[-1] == 'z':
        # compressed file names should end in 'z'
        with gzip.open(fn, mode='rt') as f:
            lines = f.readlines()
    else:
        # Someone manually uncompressed the file, probably to look at it
        with open(fn, mode='rt') as f:
            lines = f.readlines()

    curr_alignment_block = 0

    for line in lines:

        if '###' == line[0:3]:
            # We have a hit on the start of an alignment block
            curr_alignment_block += 1
            if curr_alignment_block > count:
                # We're done if we have the user's target number of rows
                break

            search_results = prog.search(line)
            alignment_header_data = list(search_results.groups())

            curr_entry['Aln_num'] = alignment_header_data[0]
            curr_entry['Prot_ID'] = alignment_header_data[1]
            curr_entry['naln'] = alignment_header_data[2]
            curr_entry['score'] = alignment_header_data[3]
            curr_entry['sid'] = alignment_header_data[4]

        elif 'Translation:' in line:
            # This is one of the header lines for the alignment block, so we can skip it
            continue
        elif 'Rotation:' in line:
            # This is one of the header lines for the alignment block, so we can skip it
            continue
        elif 'Index ' in line:
            # This is one of the header lines for the alignment block, so we can skip it
            continue
        else:
            # We have a row of actual data, so split out the data, convert the 23rd elelemnt of line to a binary value (0 if ' ' and 1 if '*'), and append that to the alignment block info
            # converting conservation symbols (* is identical restype; : is similar restype) to float values for visualization purposes
            if line[59] == '*':
                line = line[:59] + '1.0'
            elif line[59] == ':':
                line = line[:59] + '0.5'
            else:                        
                line = line[:59] + '0.0'
            
            curr_record = line.strip().split()

            curr_entry['Ind']   = curr_record[0]
            curr_entry['Ch1']   = curr_record[1]
            curr_entry['Res1']  = curr_record[2]
            curr_entry['AA1']   = curr_record[3]
            curr_entry['Ch2']   = curr_record[4]
            curr_entry['Res2']  = curr_record[5]
            curr_entry['AA2']   = curr_record[6]
            curr_entry['Dist']  = curr_record[7]
            curr_entry['Cos']   = curr_record[8]
            curr_entry['S']     = curr_record[9]

            rows.append(curr_entry.copy())

    df = pd.DataFrame(rows)

    # Converts from objects to strings
    df = df.convert_dtypes()

    # Convert numeric types from strings
    df[['Aln_num','naln','score','sid','Ind','Res1','Res2','S','Dist','Cos']] = \
        df[['Aln_num','naln','score','sid','Ind','Res1','Res2','S','Dist','Cos']].apply(pd.to_numeric)

    return df


if __name__ == '__main__':
    # test harness for these functions

    # Where the score and alignment files of interest are.
    base_path = Path('/Users/may/Projects/data/PSP/desulfovibrio_vulgaris/out/'
                     'WP_010938264.1/sadlsa_pdb70_210310')

    print('Reading score file')
    score_df = parse_sadlsa_score_file(
        str(base_path / 'WP_010938264.1_sco.dat.gz'))

    print('Reading alignment file')
    align_df = parse_sadlsa_aln_file(str(base_path / 'WP_010938264.1_aln.dat'))

    pass

