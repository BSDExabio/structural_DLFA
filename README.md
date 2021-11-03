# structural_DLFA
Deep learning and structure based hypothesis generation for functional annotation

A pipeline for bringing together hypotheses from SAdLSA and AlphaFold to uncover 
protein structure.

----

* `database` -- contains `Database`, which is a wrapper for our `sqlite3` 
  database
* `notebooks` -- jupyter notebooks that demonstrate how to use some this 
  software as well as create various visualizations; also some notebooks are 
  for working out various issues.
* `output` -- CSV output from earlier SAdLSA runs
* `parsers` -- this package contains modules for parsing Alphafold, SAdLSA, 
  and GenBank data.
* `scripts` -- scripts for ingesting data into `sqlite3` database as well as 
  lookup proteins in Genbank files.
* `seq` -- data saved from an earlier SAdLSA run
* `snippets` -- example code from which to crib
* `visualizers` -- code for visualizing Alphafold and SAdLSA data
* `website` -- support for our Missouri web site for querying and 
  visualizing `sqlite3` data.
