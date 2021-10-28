#!/usr/bin/env python3
"""
    For importing AlphaFold data into pandas dataframes.
    This is used by another utility for importing into sqlite3 database.
"""
from pathlib import Path
import pandas as pd

def parse_af_pickle_file(fns,n_models):
    """
    """
    import glob
    import pickle
    
    # grab appropriate pkl files
    pickle_files = glob.glob(fns)
    # sort pkl files
    pickle_files.sort()

    # check that not too many models are wanted than there are models available
    if n_models > len(pickle_files):
        n_models = len(pickle_files)
    
    # fill in the dictionary with model's associated plddt results
    plddt_dict = {}
    for i,model_pkl in enumerate(pickle_files[:n_models]):
        # unpickle a model's result pickle file
        with open(model_pkl,'rb') as file_:
            model_results = pickle.load(file_)
        # the pickle file's plddt dictionary element contains a numpy array of floats
        plddt_dict[i] = model_results['plddt']

    # create the dataframe
    # dataframe columns are zero-indexed numbers that can be used to index the 
    # sorted glob.glob list of files (either the pkl or model pdb files)
    df = pd.DataFrame(plddt_dict)

    return df

