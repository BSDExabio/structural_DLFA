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
