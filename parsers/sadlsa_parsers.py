#!/usr/bin/env python3
"""
    For importing SAdLSA data into pandas dataframes.

    This is used by another utility for importing into sqlite3 database.
"""

import pandas as pd


def parse_sadlsa_score_file(fn, count=10):
    """ translate given SAdLSA score file into pandas dataframe

    :param fn: string filename of score file
    :param count: how many of the top scores we want to get
    :return: pandas dataframe of
    """
    import gzip

    columns = ["seqname", "alnlen", "range1", "tmscore1", "range2", "tmscore2",
               "alnscore", "seq_id", "desc"]
    df = pd.DataFrame(columns=columns)

    with gzip.open(fn, mode='rt') as f:
        lines = f.readlines()

    for line in lines[3:count + 3]:
        temp = line.strip().split('\t|')
        rec = temp[0].split()
        rec.append(temp[1])

        series = pd.Series(rec, index=df.columns)
        df = df.append(series, ignore_index=True)

    return df


def parse_sadlsa_aln_file(fn, count=10):
    """ translate given SAdLSA alignment file into pandas dataframe

    :param fn: string filename of score file
    :param count: how many of the top scores we want to get
    :return: pandas dataframe of
    """
    import gzip
    import re
    import itertools

    prog = re.compile(
        r'### Alignment (\d+) to: (......) naln=(\d+) score=(\d+\.\d*) tms1=(\d+\.\d*) tms2=(\d+\.\d*) sid=(\d+\.\d*)')

    columns = ["Aln_num", "Prot_ID", "naln", "score", "tms1", "tms2", "sid",
               "Ind", "Res1", "AA1", "Res2", "AA2", "MeanDist", "Bin", "Prob<3",
               "Prob<5", "Prob<8"]
    df = pd.DataFrame(columns=columns)

    with gzip.open(fn) as f:
        lines = f.readlines()

    i = 0
    for line in lines:
        if i == count:
            break

        res = prog.search(line)
        if res:
            i += 1
            res = list(res.groups())
            for subline in itertools.islice(lines, int(res[2])):
                rec = res + subline.strip().replace('*', '').split()
                series = pd.Series(rec, index=df.columns)
                df = df.append(series, ignore_index=True)

    return df


if __name__ == '__main__':
# test harness for these functions
