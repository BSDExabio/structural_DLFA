#!/usr/bin/env python3
"""
    For parsing TMalign data.  It expects data like this:

/gpfs/alpine/bif135/proj-shared/rbd_work/databases/PDB70/pdb70_2022_03_19/rbd_work/pdb70_2022_03_19_structures/structure_dir/6CVL_C.pdb,1.84,201,229.0,0.88695,343.0,0.60308,0.88695
/gpfs/alpine/bif135/proj-shared/rbd_work/databases/PDB70/pdb70_2022_03_19/rbd_work/pdb70_2022_03_19_structures/structure_dir/3TIF_A.pdb,1.95,201,229.0,0.88572,230.0,0.88213,0.88572
/gpfs/alpine/bif135/proj-shared/rbd_work/databases/PDB70/pdb70_2022_03_19/rbd_work/pdb70_2022_03_19_structures/structure_dir/5LIL_B.pdb,1.96,201,229.0,0.88397,604.0,0.34851,0.88397

"""
from pathlib import Path
import pandas as pd
import gzip

def parse_tmalign_score_file(fn):
    """ Read given TMalign scores file into a pandas dataframe
    These files are written in csv format to begin with so will
    be easy to read into a pandas dataframe

    :param fn: path to the score file
    :return: pandas dataframe
    """
    
    df = pd.read_csv(fn)
    
    return df



def parse_tmalign_align_file(fn):
    """ Read given TMalign alignment file into a pandas dataframe

    :param fn: of score file
    :return: pandas dataframe
    """
    # TODO implement
    pass


if __name__ == '__main__':
    # test harness for parsers

    df = parse_tmalign_score_file('/Users/may/Projects/data/Alphafold/TMalign/WP_010939197.1/TMalign_cp_pdb70_2022_03_19/ranked_alignment_results.dat')

    pass
