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
        """
        if not date:
            date = datetime.now().date()
        if not version:
            version = '0.0'
        with self.engine.connect() as conn:
            # TODO finish creating the table for institution
            conn.execute(text('CREATE TABLE info (institution varchar(255), '
                              'date varchar(10),'
                              'version varchar(20)'))
            conn.execute(text('INSERT INTO info (institution, date, version) '
                              'VALUES (:info, :date, :version)'),
                         [{'institution' : 'Oak Ridge National Laboratory',
                           'date' : date,
                           'version' : version}])
            conn.commit()


if __name__ == '__main__':
    # test harness
