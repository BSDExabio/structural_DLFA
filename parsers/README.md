This directory contains parsers for various bioinformatics file formats.

* `af_parsers.py` -- import Alphafold data into a dataframe
* `apoc_parsers.py` -- import APOC data into a dataframe
* `genbank.py` -- this will read a Genbank file and return mappings between 
  protein IDs and their corresponding locus tags; note that locus tags will 
  have their underscores stripped
* RCSB data
  * `rcsb_query.py` -- query RCSB site for Uniprot ID for a given protein and chain ID
  * `rcsb_scraper.py` -- similar to `rcsb_query.py`, but scrape the web site rendered pages, instead
      > Prefer `rcsb_query.py` since that puts less of a load on RCSB servers
* `sadlsa_parser.py` -- for reading SAdLSA files
* `tmalign_parser.py` -- parse TMalign-APOC data into dataframes
