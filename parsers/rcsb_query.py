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
    uniprot_id = None

    r = requests.post('https://data.rcsb.org/graphql',
                      json={'query' : QUERY_STR,
                            'variables' : {'id' : protein.upper()}})
    data = r.json()

    if 'error' in data:
        # There was something wrong with the query
        # TODO extract some meaningful error from this for an informative
        #  exception
        return None

    chain = chain.upper() # Normalize to chain IDs being upper case

    for entry in data['data']['entries']:
        for entity in entry['polymer_entities']:
            if chain in entity['rcsb_polymer_entity_container_identifiers']['asym_ids'] or \
                chain in entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']:
                uniprot_id = entity['uniprots'][0]['rcsb_id']
                break

    return uniprot_id


if __name__ == '__main__':
    results = query_uniprot('6lzm', 'A')
    assert results == 'P00720'


    pass
