import pandas as pd

def parse_sadlsa_score_file(fn, count=10):
    import gzip

    columns = ["seqname", "alnlen", "range1", "tmscore1", "range2", "tmscore2", "alnscore", "seq_id", "desc"]
    df = pd.DataFrame(columns = cols)

    with gzip.open(fn) as f:
        lines=f.readlines()

    for line in lines[3:count+3]:
        temp = line.strip().split('\t|')
        rec = temp[0].split()
        rec.append(temp[1])

        series = pd.Series(rec, index=df.columns)
        df = df.append(series, ignore_index=True)

    return df

def parse_sadlsa_aln_file(fn, count=10):
    import gzip
    import re
    import itertools

    prog = re.compile(r'### Alignment (\d+) to: (......) naln=(\d+) score=(\d+\.\d*) tms1=(\d+\.\d*) tms2=(\d+\.\d*) sid=(\d+\.\d*)')

    columns = ["Aln_num", "Prot_ID", "naln", "score", "tms1", "tms2", "sid", "Ind", "Res1", "AA1", "Res2", "AA2", "MeanDist", "Bin", "Prob<3", "Prob<5", "Prob<8"]
    df = pd.DataFrame(columns = columns)


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
