#!/usr/bin/env bash
#
# Will create sqlite3 database from a CSV file that was read as a
# pandas dataframe, and then do an external query to show that it worked.

# create 'save_pandas.db' sqlite3 database
python3 ./pandas_2_sqlite3.py

# do a stupid query on it to show that all the data is there  :P
