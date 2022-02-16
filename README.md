# structural_DLFA
Deep learning and structure based hypothesis generation for functional annotation

A pipeline for bringing together hypotheses from SAdLSA and AlphaFold to uncover 
protein structure.

----

## Directories

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

## Files

* `requirements.txt` -- python dependencies for `database`, `scripts`, and 
  `parsers`
* `Dockerfile` -- for building web site container
* `docker-compose.yml` -- for spinning up the container

----

## Spinning up web site

1. Ensure your environment is setup
   1. docker is installed
      > We recommned installing via the website and not via home brew or a 
      linux package manager.
   2. git lfs is installed
   3. you have checked out `website/static/dlfa.db` via `git lfs`
2. `docker-compose up`
   > You will see logging messages in your terminal for all GET and POST 
   > transactions
3. Surf to `https://0.0.0.0:5000`

(Yes, eventually the port and host IP will change.)
