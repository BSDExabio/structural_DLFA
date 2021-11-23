#!/usr/bin/env python3
"""
    Main support for the website
"""
import sys
from pathlib import Path
from markupsafe import escape
from database.database import Database

from flask import Flask, render_template, request
app = Flask(__name__)

sqlite3_db = Path('../db/dlfa.db')
db = Database(sqlite3_db)

@app.route('/')
def index():
    info = db.get_info_table()
    score_proteins = db.query_unique_sadlsa_score_proteins()
    alignment_proteins = db.query_unique_sadlsa_alignments_proteins()
    return render_template('index.html', info=info,
                           num_score_proteins=len(score_proteins),
                           num_alignment_proteins=len(alignment_proteins))

@app.route('/align_output', methods=['POST'])
def align_output():
    protein = request.form.get('protein')
    results = db.query_sadlsa_alignments(protein)
    return render_template('align_output.html', results=results)

@app.route('/user/<username>')
def show_user_profile(username):
    return f'User {escape(username)}'

@app.route('/post/<int:post_id>')
def show_post(post_id):
    return f'Post {post_id}'

@app.route('/path/<path:subpath>')
def show_subpath(subpath):
    return f'Subpath {escape(subpath)}'

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
