#!/usr/bin/env python3
"""
    Query RCSB site Uniprot ID associated with a given protein and chain ID.

    (Written by Mark Coletti, colettima@ornl.gov.)
"""
import requests

QUERY_STR = \
'''
query($id: String!)
{
  entries(entry_ids:[$id]){
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
      }
      uniprots {
        rcsb_id
      }
    }
  }
}
'''

def query_uniprot(protein, chain):
    """ Get Uniprot ID for the given protein and protein chain

    :param protein: for which to get Uniprot ID
    :param chain: specific chain in that protein
    :return: Uniprot ID or None if not found
    """
    r = requests.post('https://data.rcsb.org/graphql',
                      json={'query' : "query($id: String!){entry(entry_id:$id){exptl{method}}}",
                            'variables' : {'id' : protein.upper()}})
    pass


if __name__ == '__main__':
    results = query_uniprot('6lzm', 'A')

    pass