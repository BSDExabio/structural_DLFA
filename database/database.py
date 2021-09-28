#!/usr/bin/env python3
"""
    Database class for managing and accessing sqlite3 database.
"""
from datetime import datetime
from sqlalchemy import create_engine, text
import json
import sqlite3

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


    def query_enzyme(self, class_num, subclass_num, subsubclass_num):
        """ Queries the enzyme table

        We return a JSON object because that works well with web sites.

        :param class_num: class number of enzyme
        :param subclass_num: subclass number of enzyme
        :param subsubclass_num: subsubclass number of enzyme

        :return: JSON of enzyme or None if not found
        """
        with self.engine.connect() as conn:
            result = conn.execute(text("SELECT class_name, "
                                       "subclass_desc, "
                                       "subsubclass_desc "
                                       "FROM enzyme_desciptions "
                                       "WHERE class = :class_num AND "
                                       "subclass = :subclass_num AND "
                                       "subsubclass = :subsubclass_num" ),
                                  {'class_num' : class_num,
                                   'subclass_num' : subclass_num,
                                   'subsubclass_num' : subsubclass_num})

            for class_name, subclass_desc, subsubclass_desc in result:
                return json.dumps({'class_name' : class_name,
                                   'subclass_desc' : subclass_desc,
                                   'subsubclass_desc' : subsubclass_desc})

            # If we got this far, then there were no hits.
            return None


if __name__ == '__main__':
    # test harness
    dlfa_db = Database('/tmp/dlfa.db')
    dlfa_db.add_info_table()

    try:
        enzyme = dlfa_db.query_enzyme(1, 1, 1)
        print(f'Got {enzyme}')
    except sqlite3.OperationalError as e:
        print(f'{e}')
