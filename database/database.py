#!/usr/bin/env python3
"""
    Database class for managing and accessing sqlite3 database.
"""
from datetime import datetime
from sqlalchemy import create_engine, text
import json
import sqlite3


def write_df_to_db(data_frame, table, database):
    """ write the given EC dataframe to the sqlite3 db

    This will *append* data not do an upsert, so if you run this multiple
    times you with the same protein, it will appear more than once in the
    database.

    TODO there has to be a way to do an upsert with sqlachemy and sqlite3.

    :param data_frame: dataframe of EC nomenclature
    :param table: to write to in database
    :param database: name of sqlite3 database
    :return: none.
    """
    # this will create the sqlite3 database if it does not already exist;
    # echo=True to have it be noisy about what it's doing; turn that off for
    # production, natch
    engine = create_engine(f'sqlite:///{database}', echo=False)

    with engine.connect() as sqlite_connection:
        # if_exists='replace' means dropping the whole table first and then
        # replacing if with new values.  See
        # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas
        # .DataFrame.to_sql.html for more information.  (You may want to
        # append, instead.  Also, I was hoping there would be support for
        # "upsert", which I don't see; i.e., add new information if not
        # already there, but replace it if it is.)
        data_frame.to_sql(table, sqlite_connection,
                          index=False,
                          if_exists='append')


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

        >>> dlfa_db = Database('/tmp/dlfa.db')
        >>> enzyme = dlfa_db.query_enzyme(1, 1, 1)
        >>> print(enzyme)
        {"class_name": "Oxidoreductases", "subclass_desc": "Acting on the CH-OH group of donors", "subsubclass_desc": "With NAD(+) or NADP(+) as acceptor"}

        :param class_num: class number of enzyme
        :param subclass_num: subclass number of enzyme
        :param subsubclass_num: subsubclass number of enzyme

        :return: JSON of enzyme or None if not found
        """
        with self.engine.connect() as conn:
            rows = conn.execute(text("SELECT class_name, "
                                       "subclass_desc, "
                                       "subsubclass_desc "
                                       "FROM enzyme_desciptions "
                                       "WHERE class = :class_num AND "
                                       "subclass = :subclass_num AND "
                                       "subsubclass = :subsubclass_num" ),
                                  {'class_num' : class_num,
                                   'subclass_num' : subclass_num,
                                   'subsubclass_num' : subsubclass_num}).fetchall()

            if rows != []:
                # We should only get a single row for this kind of query
                class_name, subclass_desc, subsubclass_desc = rows[0]
                return json.dumps({'class_name' : class_name,
                                   'subclass_desc' : subclass_desc,
                                   'subsubclass_desc' : subsubclass_desc})

            # If we got this far, then there were no hits.
            return None

    def query(self, sql):
        """ Do SQL query on database

        >>> dlfa_db = Database('/tmp/dlfa.db')
        >>> results = dlfa_db.query('select class_name, subsubclass_desc from enzyme_desciptions where class = 1 and subclass = 1 and subsubclass = 1;')
        >>> print(results)
        [["Oxidoreductases", "With NAD(+) or NADP(+) as acceptor"]]

        :param sql: SQL query
        :returns: JSON formatted query results
        """
        with self.engine.connect() as conn:
            rows = conn.execute(text(sql))

            if rows == []:
                return None
            else:
                return json.dumps([tuple(row) for row in rows])



if __name__ == '__main__':
    # test harness
    dlfa_db = Database('/tmp/dlfa.db')
    dlfa_db.add_info_table()

    try:
        enzyme = dlfa_db.query_enzyme(1, 1, 1)
        print(f'Got {enzyme}')
        enzyme = dlfa_db.query_enzyme(99, 99, 99)
        print(f'Got {enzyme}')
    except sqlite3.OperationalError as e:
        print(f'{e}')
