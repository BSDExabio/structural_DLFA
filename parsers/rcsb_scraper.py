#!/usr/bin/env python3
"""
    Scrape RCSB site for chains and enzymes associated with a given protein.
"""
import re
from pprint import pprint

import requests
from bs4 import BeautifulSoup


def scrape_rcsb(protein):
    """ Get chains and enzymes for given protein

    :param protein: for which to get chains and enzymes
    :return: chains and their enzymes, or None if none found
    """
    def extract_chain(table_row):
        chain = table_row.contents[1].a.text

        enzyme_content = table_row.contents[4].find('a', class_='querySearchLink',
                                          href=re.compile('rcsb_ec_lineage'))

        if enzyme_content:
            # FIXME, need to accommodate multiple enzymes
            enzyme = enzyme_content.text
        else:
            return None

        return chain, enzyme


    scrape_results = requests.get(f'https://www.rcsb.org/structure/{protein.upper()}')
    if scrape_results.status_code != 200:
        return None

    result = {'protein' : protein, 'chains' : []}

    soup = BeautifulSoup(scrape_results.text, 'html.parser')

    # Find all the relevant table rows
    table_rows = soup.find_all(id=re.compile('macromolecule-entityId'))

    for table_row in table_rows:
        new_chain = extract_chain(table_row)
        if new_chain: # found a chain with one or more enzymes
            result['chains'].append(new_chain)

    return result



if __name__ == '__main__':
    # results = scrape_rcsb('6lzm')
    # pprint(results)
    results = scrape_rcsb('5fvk')
    pprint(results)
