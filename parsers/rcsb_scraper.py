#!/usr/bin/env python3
"""
    Scrape RCSB site for chains and enzymes associated with a given protein.
"""
import requests

def scrape_rcsb(protein):
    """ Get chains and enzymes for given protein

    Returned dict will have this structure:

    {'protein' : {'chain1' : {'enzyme1' : 'source', 'enzyme2' : 'source'},
                  'chain2' : {'enzyme3' : 'source'}}}

    You could stack these in a list to later build a pandas dataframe from.

    :param protein: for which to get chains and enzymes
    :return: dict of chains and enzymes, or None if an error
    """
    result = requests.get(f'https://www.rcsb.org/structure/{protein.upper()}')
    if result.status_code != 200:
        return None

if __name__ == '__main__':
    results = scrape_rcsb('6lzm')

    pass
