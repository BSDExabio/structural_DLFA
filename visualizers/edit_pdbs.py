#!/usr/bin/env python3

"""
"""

def grab_structure(pdbID):
    """
    Pull the structure from the RCSB database.
    INPUT:
    :param pdbID: string; 4 character long code associated with a structure on the RCSB
    OUTPUT:
    :return: a mmtf structure object
    """

    import mmtf

    return mmtf.fetch(pdbID)


def make_atm_sel(universe,sel_string):
    """
    create a MDAnalysis atom group based on the given selection string.
    INPUT:
    :param universe: MDAnalysis Universe object; universe that will be used to create the atom group
    :param sel_string: string; a MDAnalysis atom selection string used to create the atom group
    OUTPUT:
    :return: a MDAnalysis atom group object
    """
    return universe.select_atoms(sel_string)


def create_2d_array(array):
    """
    create a 2d array from a 1d array, where the first column in the new 2d array is a zero-indexed count of the rows
    used to create a mirror of SAdLSA df 2d arrays that are fed into edit_pdbs.edit_pdb function
    INPUT:
    :param array: 1d array of N size
    OUTPUT:
    :return: a 2d numpy array of (N,2) size with [:,0] as zero-indices and [:,1] as the array values
    """
    import numpy as np
    return np.c_[np.arange(array.shape[0]),array]


def make_zero_indexed(array):
    """
    subtract one from every element in array; used to change the array from one-indexing to zero-indexing to make uniform with python indexing.
    used to fix residue indices in the SAdLSA arrays so they can be used to index residues in MDAnalysis Universes
    INPUT:
    :param array: 1d array
    OUTPUT:
    :return: 1d array with one subtracted from every element
    """
    return array - 1


def edit_pdb(target_structure, metric_data, out_filename, default_value = -1.00, working_dir = './', sel_string = 'protein'):
    """
    Editing the atom data of a pdb-formatted data set; filling the 61-66 columns of the pdb to the float values of the metric_data.

    INPUTS:
    :param target_structure: file or object that can be used to initiate a MDAnalysis universe object; current use cases: mmtf structure object or string that points to a .pdb file.
    :param metric_data: a 2d array the contains the metric of interest; first column corresponds to the one-indexed residue ids; 2nd column corresponds to the float values associated with the metric of interest. 
    :param out_filename: a string that acts as a descriptor for the metric; should not contain any spaces.
    :param default_value: float; this value is used in b-factor columns for residues that do not have a set value. Default: -1.00. 
    :param working_dir: string; local or global path string where files will be saved; path needs to exist already. Default: './' 
    :param sel_string: string; MDAnalysis atom selection string used to create the important group that will have b-factors changed. Default: 'protein' 
    :return: a string associated with the path to the new saved pdb file; if no file is saved, returns 0 (int).

    OUTPUTS: 
    A new .pdb file that contains the metric of interest in the b-factor column. If the system is too large, no file is written. 
    """

    import MDAnalysis
    import numpy as np
    import warnings

    # setting output file name
    file_name = working_dir + out_filename
    # loading the mmtf structure object into an MDAnalysis universe object
    u = MDAnalysis.Universe(target_structure)
    # creating the atom selection groups for the structure
    sel = make_atm_sel(u,sel_string)
    #_all = make_atm_sel(u,'all')

    # moving all atoms of chain of interest to have the CoM to be the origin
    #_all.translate(-sel.center_of_mass())
    sel.translate(-sel.center_of_mass())
    
    # creating the one-indexed array of sel's residue indices of resolved 
    # residues; used to mirror the one-indexed SAdLSA resid column
    nRes_range = range(0,sel.n_residues)
    for i in nRes_range:
        if i in metric_data[:,0]:
            idx = np.argwhere(i == metric_data[:,0])
            sel.residues[i].atoms.tempfactors = float(metric_data[idx,1])   # will be printed with 2 decimal points; 
        else:
            sel.residues[i].atoms.tempfactors = default_value
    
    with warnings.catch_warnings():
        # ignore some annoying warnings from sel.write line due to missing information (chainIDs, elements, and record_types). 
        warnings.simplefilter('ignore',UserWarning)
        if sel.n_atoms > 99999.:
            print('Number of atoms is too large for pdb file format; need to visualize these results without writing to file.')
            return 0
        else:
            #_all.write(file_name)
            sel.write(file_name)
            return file_name


def simple_align(pdb_list,aln_selection = 'protein and name CA'):
    """
    use the aln_selection atom group to align files in the pdb_list; overwrites the aligned structures to their respective file name
    :param pdb_list: list of strings; each element is a global or local path to a .pdb file
    :param aln_selection: string; a MDAnalysis atom selection string; default: 'protein and name CA'
    OUTPUTS:
    Overwrites the provided .pdb files with the aligned coordinates; maintains all other information
    """
    import MDAnalysis
    from MDAnalysis.analysis.align import rotation_matrix
    import warnings
    
    universe = []
    all_sels = []
    aln_sels = []
    with warnings.catch_warnings():
        # ignore some annoying warnings from sel.write line due to missing information (chainIDs, elements, and record_types). 
        warnings.simplefilter('ignore',UserWarning)
        for pdb in pdb_list:
            temp = MDAnalysis.Universe(pdb)
            temp_all = temp.select_atoms('all')
            temp_aln = temp.select_atoms(aln_selection)
            # remove translational differences in structures
            temp_all.translate(-temp_aln.center_of_mass())
            all_sels.append(temp_all)
            aln_sels.append(temp_aln)

        pos0 = aln_sels[0].positions
        for i in range(len(pdb_list)):
            # calculate the rotation matrix needed to minimize the RMSD between alignment selections of the two universe
            R, rmsd = rotation_matrix(aln_sels[i].positions,pos0)
            # apply the matrix to remove rotational differences in the two structures
            all_sels[i].rotate(R)
            all_sels[i].write(pdb_list[i])


