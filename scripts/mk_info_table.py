#!/usr/bin/env python3
"""
    Create the sqlite3 info table that describes the version and date
    for the database.
"""
import argparse
from database.database import Database

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='For creating or updating '
                                                 'the sqlite info table')
    parser.add_argument('database', help='Database to be updated')
    parser.add_argument('--date', help='Date of this version; defaults to today')
    parser.add_argument('--version', help='Version of the database; defaults to 0.0')

    args = parser.parse_args()

    dlfa_db = Database(args.database)
    dlfa_db.add_info_table(version=args.version, date=args.date)
