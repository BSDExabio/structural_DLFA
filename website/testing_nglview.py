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

pdb_file = '2HW4_A_S.pdb'

#@app.route('/')
#def index():
#    info = db.get_info_table()
#    score_proteins = db.query_unique_sadlsa_score_proteins()
#    alignment_proteins = db.query_unique_sadlsa_alignments_proteins()
#    return render_template('index.html', info=info,
#                           num_score_proteins=len(score_proteins),
#                           num_alignment_proteins=len(alignment_proteins))
#

#@app.route('/')
#def tdmol_visualize():
#    import nglview as nv
#    print("nglview version = {}".format(nv.__version__))
#
#    #u = MDAnalysis.Universe(pdb_file)
#    #sel = u.select_atoms('protein')
#    #w = nv.show_mdanalysis(sel)
#    #w.clear_representations()
#    #w.add_representation('cartoon',color='bfactor')
#    #w.camera = 'orthographic'
#    return render_template("temp.html", pdb = pdb_file.split('/')[-1][:4], pdb_file = pdb_file)    # view=w, 

@app.route('/')
def nglv_visualize():
    #import nglview as nv
    #print("nglview version = {}".format(nv.__version__))

    #u = MDAnalysis.Universe(pdb_file)
    #sel = u.select_atoms('protein')
    #w = nv.show_mdanalysis(sel)
    #w.clear_representations()
    #w.add_representation('cartoon',color='bfactor')
    #w.camera = 'orthographic'
    return render_template("align_visualization.html", align_type = 'SAdLSA', seqid = 'Test', pdb = pdb_file.split('/')[-1][:4], metric = 'S', pdb_file = pdb_file)    # view=w, 



if __name__ == '__main__':
    #if not sqlite3_db.exists():
    #    print(f'{sqlite3_db} does not exist; please create and restart')
    #    sys.exit(1)

    app.run(debug=True)
