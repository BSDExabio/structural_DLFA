#!/usr/bin/env python3
"""
    For importing AlphaFold data into pandas dataframes.
    This is used by another utility for importing into sqlite3 database.
"""
import pandas as pd

def parse_af_pickle_files(fn_list):
    """ translate given AlphaFold result files into a pandas dataframe
    saving the pLDDT results for all user-specified models
    :param fn_list: list of strings; the global or local paths to pickle files
    :return: pandas dataframe of models' pLDDT values
    """
    import pickle

    # fill in the dictionary with model's associated plddt results
    plddt_dict = {}
    for i, model_pkl in enumerate(fn_list):
        # unpickle a model's result pickle file
        with open(model_pkl,'rb') as file_:
            model_results = pickle.load(file_)
        # the pickle file's plddt dictionary element contains a numpy array of floats
        plddt_dict[i] = model_results['plddt']

    # create the dataframe
    # dataframe columns are zero-indexed numbers that can be used to index the 
    # list of files (either the pkl or model pdb files)
    df = pd.DataFrame(plddt_dict)

    return df

