#!/usr/bin/env python3
"""
    Scrape RCSB site for chains and enzymes associated with a given protein.

    (Written by Mark Coletti, colettima@ornl.gov.)
"""
import re
from pprint import pprint

import requests
from bs4 import BeautifulSoup


def scrape_rcsb(protein, chain):
    """ Get enzymes for given protein and chain

    :param protein: for which we are searching
    :param chain: specific chain of that protein
    :return: list of enzymes, or None if there is no match
    """
    def extract_chain_enzymes(table_row, chain):
        # Snip out all the chain IDs in the first table row data cell.  They
        # will all have HTML 'a' anchor tags.
        found_chains = [x.text for x in table_row.contents[1].find_all('a')]

        if chain not in found_chains:
            # This does not contain the chain we're looking for
            return None

        # Now snip out the EC, which will be an HTML anchor with the
        # querySearchLink class.
        enzyme_content = table_row.contents[4].find('a',
                                                    class_='querySearchLink',
                                          href=re.compile('rcsb_ec_lineage'))

        if enzyme_content:
            # FIXME, need to accommodate multiple enzymes
            enzyme = enzyme_content.text
        else:
            return None

        return enzyme


    scrape_results = requests.get(f'https://www.rcsb.org/structure/{protein.upper()}')
    if scrape_results.status_code != 200:
        return None

    result = [] # will contain zero or more enzymes

    soup = BeautifulSoup(scrape_results.text, 'html.parser')

    # Find all the relevant table rows
    table_rows = soup.find_all(id=re.compile('macromolecule-entityId'))

    for table_row in table_rows:
        new_enzymes = extract_chain_enzymes(table_row, chain)
        if new_enzymes: # found one or more enzymes
            result.append(new_enzymes)

    return result



if __name__ == '__main__':
    # Simple case of one chain to one enzyme
    results = scrape_rcsb('6lzm', 'A')
    pprint(f'6lzm, chain A: {results}')
    assert results == ['3.2.1.17']

    results = scrape_rcsb('5fvk', 'B')
    pprint(f'5fvk, chain B: {results}')
    assert results == ['3.6.4.6']
