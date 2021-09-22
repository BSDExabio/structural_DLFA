#!/usr/bin/env python3
"""
    This demonstrates how to insert a new sqlite3 table from a pandas
    dataframe.

    Cribbed from:

    https://www.fullstackpython.com/blog/export-pandas-dataframes-sqlite-sqlalchemy.html
"""
import pandas as pd
from sqlalchemy import create_engine

if __name__ == '__main__':

    # We're just going to arbitrarily read and insert ../output/sadlsa_alignment_quality.csv
    # as a sqlite3 table.
    sadlsa_align_qual_df = pd.read_csv('../output/sadlsa_alignment_quality.csv')

    # this will create the sqlite3 database if it does not already exist;
    # echo=True to have it be noisy about what it's doing; turn that off for
    # production, natch
    engine = create_engine('sqlite:///save_pandas.db', echo=True)

    with engine.connect() as sqlite_connection:
        # What we wanna call our table
        sqlite_table = 'sadlsa_alignment_quality'

        # if_exists='replace' means dropping the whole table first and then
        # replacing if with new values.  See
        # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas
        # .DataFrame.to_sql.html for more information.  (You may want to
        # append, instead.  Also, I was hoping there would be support for
        # "upsert", which I don't see; i.e., add new information if not
        # already there, but replace it if it is.)
        sadlsa_align_qual_df.to_sql(sqlite_table, sqlite_connection,
                                    if_exists='replace')

    pass
