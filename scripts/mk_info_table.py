#!/usr/bin/env python3
"""
    Create the sqlite3 info table that describes the version and date
    for the database.

    usage: mk_info_table.py [-h] [--date DATE] [--version VERSION] database

    For creating or updating the sqlite info table

    positional arguments:
      database           Database to be updated

    optional arguments:
      -h, --help         show this help message and exit
      --date DATE        Date of this version; defaults to today
      --version VERSION  Version of the database; defaults to 0.0
"""
import argparse
from database.database import Database

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='For creating or updating '
                                                 'the sqlite info table')
    parser.add_argument('database', help='Database to be updated')
    parser.add_argument('--date', help='Date of this version; defaults to '
                                       'today')
    parser.add_argument('--version', help='Version of the database; defaults '
                                          'to 0.0')

    args = parser.parse_args()

    dlfa_db = Database(args.database)
    dlfa_db.add_info_table(version=args.version, date=args.date)
