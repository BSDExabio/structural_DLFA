#!/usr/bin/env python3
"""
    Database class for managing and accessing sqlite3 database.
"""
from datetime import datetime
from sqlalchemy import create_engine, text

class Database():
    """
        For managing and accessing the DLFA sqlite3 database.
    """
    def __init__(self, db_filename):
        """
        :param db_filename: file for sqlite3 database
        """
        self.db_filename = db_filename

        self.engine = create_engine(f'sqlite:///{db_filename}')


    def add_info_table(self, version=None, date=None):
        """
            This will add a table in the sqlite3 database that will identify
            the source of this database.

            TODO need to drop table before hand if run on existing database
            so tht we only have the one canonical row of info.
        """
        if not date:
            date = datetime.now().date().isoformat()
        if not version:
            version = '0.0'
        with self.engine.connect() as conn:
            conn.execute(text('CREATE TABLE IF NOT EXISTS '
                              'info (institution TEXT, '
                              'date TEXT, '
                              'version TEXT)'))
            conn.execute(text('INSERT INTO info (institution, date, version) '
                              'VALUES (:institution, :date, :version)'),
                         [{'institution' : 'Oak Ridge National Laboratory',
                           'date' : date,
                           'version' : version}])


if __name__ == '__main__':
    # test harness
    dlfa_db = Database('/tmp/dlfa.db')
    dlfa_db.add_info_table()
