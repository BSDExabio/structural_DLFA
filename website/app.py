#!/usr/bin/env python3
"""
    Main support for the website
"""
import sys
from pathlib import Path

from database.database import Database

from flask import Flask, render_template, request
app = Flask(__name__)

sqlite3_db = Path('../db/dlfa.db')

# @app.route('/', methods = ['POST', 'GET'])
# def query():
#     return render_template('query.html')
#
# @app.route('/list', methods = ['POST'])
# def list():
#     protein_id = request.form['protein']
#
#     rows = Database.query_sadlsa_alignments(protein_id)



if __name__ == '__main__':
    if not sqlite3_db.exists():
        print(f'{sqlite3_db} does not exist; please create and restart')
        sys.exit(1)

    app.run(debug=True)
