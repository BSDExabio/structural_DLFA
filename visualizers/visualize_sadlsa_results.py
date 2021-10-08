#!/usr/bin/env python3
"""
    For pulling information out of SAdLSA pandas dataframes and visualizing results.

    Expected functions to be made:
    1) query_dataframe. Pull important information from the dataframe(s) created by parsing SAdLSA alignment results.
    2) gather_alignment_structure. Take a PDB ID from the pulled info and search in a structure database (local) or download the relevant structure file from RSCB (remote). Return meta and atom data for later use. 
    3) gather_alignment_data. Take 

"""
#from pathlib import Path
#import pandas as pd

def query_alignment_dataframe():
    """
    Code to pull information from a sqlite database or pandas dataframe or loading a cvs file into memory and reading from that... how should we go about doing this?
    """


def grab_structure(pdbid,working_dir='./'):
    """
    Pull the structure from the RSCB database.

    :param subclass_num: subclass number of enzyme
    :param subsubclass_num: subsubclass number of enzyme
    :return: JSON of enzyme or None if not found
    

    """

    from urllib import request

    pdb_urls = 'https://files.rcsb.org/view/%s.pdb' # static url for the pdb file associated with pdb ids; %s is replaced with the pdb id.
    fasta_urls = 'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed' # static url for the fasta file associated with pdb ids; %s is replaced with the pdb id.
    cif_urls = 'https://files.rcsb.org/view/%s.cif' # static url for the mmcif file assocaited with pdb ids; %s is replaced with the pdb id.

    urls = [pdb_urls,fasta_urls,cif_urls]
    for url in urls:
        try:
            # grab the url object associated with the url string 
            response = request.urlopen(url%(pdbid))
        except:
            print(url%(pdbid) + ' returns an error. No ' + extension + ' written out.')
            return
        # convert the url object to a string for reading
        response_str = str(response.read())
        if 'Error Page' in response_str:
            print(url%(pdbid) + ' returns an error. No ' + extension + ' written out.')
            return
        # split the long string into a list of strings associated with each line of the file
        response_list = response_str[2:-1].split('\\n')
        # parse the file's lines
        #currently only saving the file to storage.
        # in the pdb or mmcif file, I should instead grab specific REMARK lines and all the ATOM lines associated with the protein chain of interest, store them into a data structure for use later.
        with open(working_dir+pdb_code+'.'+extension,'w') as W:
            for line in response_list:
                W.write(line+'\n')








