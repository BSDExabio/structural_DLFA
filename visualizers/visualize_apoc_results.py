#!/usr/bin/env python3

"""
    Functions for pulling information out of APOC pandas dataframes and visualizing results.
"""


def query_score_df(score_df):
    """
    Code to pull information from a pandas dataframe created from the apoc score file
    :param score_df: pandas dataframe object; filled with the score results
    :return: an array of lists; each element in the array is a list containing the pdb_chainID string ("seqname" column), tmscore ("tmscore"; float), tmscore's rmsd ("tmrmsd"; float), alignment length ("tmalnlen"; float), and sequence identity ("seq_id"; float) data.
    """
    grab_columns = ["seqname", "tmscore", "tmrmsd", "tmalnlen", "seq_id"]
    return score_df.filter(grab_columns).values


def query_alignment_df(align_df, seqname, metric_string):
    """
    Code to pull information from a pandas dataframe created from the apoc alignment file
    :param align_df: pandas dataframe object; filled with the alignment results
    :param seqname: string; associated with a specific target protein structure, specified in the "seqname" column of the score results. 
    :return: a 2d numpy array; 0th column corresponds to the target protein structure's residue index, 1st column corresponds to the metric of interest 
    """
    temp_df = align_df[align_df['Prot_ID'].str.contains(seqname)]
    grab_columns = ["Res2",metric_string]
    return temp_df.filter(grab_columns).values


