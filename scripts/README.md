# Deep Learning-based Functional Annotation (DLFA) related scripts.

This directory contains general purpose scripts.

* `ingest_sadlsa.sh` -- this will ingest all SAdLSA output to an sqlite3 
  database; calls `sadlsa_2_sqlige3.py`
* `ec_2_sqlite3.py` -- import the enzyme database into sqlite3 database
* `sadlsa_2_sqlige3.py` -- import a specific protein SAdLSA run into the 
  sqlite3 database
* `lookup_protein.py` -- reads a CSV or space delimited file and will add 
  columns for the protein IDs for any specified column of locus tags
