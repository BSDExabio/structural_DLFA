#!/usr/bin/env python3

"""
    For pulling information out of SAdLSA pandas dataframes and visualizing results.

    Expected functions to be made:
    1) query_dataframe. Pull important information from the dataframe(s) created by parsing SAdLSA alignment results. (collecting the beta-column values; two dimensional array of nRes_aligned x 2 shape; first column is the resid that had successfully been aligned, second column is the metric of interest)
    2) gather_alignment_structure. Take a PDB ID from the pulled info and search in a structure database (local) or download the relevant structure file from RSCB (remote). Return meta and atom data for later use. 
    3) load pdb files and edit the beta column to have the desired values...
    4) visualize in a molecular visualizer code... I have code to save files for VMD but I don't think this is what we want...
"""


def query_score_df(score_df):
    """
    Code to pull information from a pandas dataframe created from the sadlsa score file
    :param score_df: pandas dataframe object; filled with the score results
    :return: an array of lists; each element in the array is a list containing the pdb_chainID string ("seqname" column), predicted tmscore ("tmscore1"; float), alignment score ("alnscore"; float), alignment length ("alnlen"; float), and sequence identity ("seq_id"; float) data.
    """
    grab_columns = ["seqname","tmscore1","alnscore","alnlen","seq_id"]
    return score_df.filter(grab_columns).values


def query_alignment_df(align_df, seqname, metric_string):
    """
    Code to pull information from a pandas dataframe created from the sadlsa alignment file
    :param align_df: pandas dataframe object; filled with the alignment results
    :param seqname: string; associated with a specific target protein structure, specified in the "seqname" column of the score results. 
    :return: a 2d numpy array; 0th column corresponds to the target protein structure's residue index, 1st column corresponds to the metric of interest 
    """
    temp_df = align_df[align_df['Prot_ID'].str.contains(seqname)]
    grab_columns = ["Res2",metric_string]
    return temp_df.filter(grab_columns).values


def grab_structure(pdbID):
    """
    Pull the structure from the RCSB database.
    
    INPUT:
    :param pdbID: string; 4 character long code associated with a structure on the RCSB
    :return: a path string pointing to the saved structural file
    """

    import mmtf

    return mmtf.fetch(pdbID)


def edit_pdb(target_structure, pdbchainID, metric_data, metric_type, default_value = -1.00, working_dir = './'):
    """
    Editing the atom data of a pdb-formatted data set; filling the 61-66 columns of the pdb to the float values of the metric_data.

    INPUTS:
    :param target_structure: mmtf structure object.
    :param pdbchainID: string that denotes the pdb and chain ID of the protein chain used as the alignment target.
    :param metric_data: a 2d array the contains the metric of interest; first column corresponds to the one-indexed residue ids; 2nd column corresponds to the float values associated with the metric of interest. 
    :param metric_type: a string that acts as a descriptor for the metric; should not contain any spaces.
    :return: a string associated with the path to the new saved pdb file; if no file is saved, returns 0 (int).

    OUTPUTS: 
    A new .pdb file that contains the metric of interest in the b-factor column. If the system is too large, no file is written. 
    """

    import MDAnalysis
    import numpy as np
    import warnings

    # setting output file name
    file_name = working_dir + pdbchainID + '_' + metric_type + '.pdb'
    chainID = pdbchainID.split('_')[1]
    # loading the mmtf structure object into an MDAnalysis universe object
    u = MDAnalysis.Universe(target_structure)
    # creating the atom selection groups for the structure
    prot = u.select_atoms('protein and segid %s'%(chainID))
    _all = u.select_atoms('segid %s'%(chainID))
    # moving all atoms of chain of interest to have the CoM to be the origin
    _all.translate(-prot.center_of_mass())
    
    # creating the one-indexed array of protein residue indices of resolved 
    # residues; used to mirror the one-indexed SAdLSA resid column
    nRes_range = range(1,prot.n_residues+1)
    for i in nRes_range:
        if i in metric_data[:,0]:
            idx = np.argwhere(i == metric_data[:,0])
            prot.residues[i-1].atoms.tempfactors = float(metric_data[idx,1])   # will be printed with 2 decimal points; 
        else:
            prot.residues[i-1].atoms.tempfactors = default_value
    if prot.n_atoms > 99999.:
        print('Number of atoms is too large for pdb file format; need to visualize these results without writing to file.')
        return 0
    else:
        prot.write(file_name)
        return file_name


def create_vmd_vis_state(vis_state_file_name, colorbar_file_name, pdb_file_name, metric_max, max_color, metric_min=0., min_color='lightgray', under_color='white',colorbar_label='Metric'):
    """
    Writing the VMD visualization state file to view the protein structure with the metric values in the b-factor column of the pdb
    Also, creating a colorbar figure that mirrors the color range used in the VMD vis state. 

    INPUTS:
    :param vis_state_file_name: string of a global or local path that the vis state file will be saved as
    :param colorbar_file_name: string of a global or local path that the colorbar will be saved as
    :param colorbar_file_name: string of a global or local path of the pdb to be visualized
    :param metric_max: float that sets the maximum value of the color scale
    :param max_color: string of a color to be used to label the maximum value; must be a color accepted by matplotlib.colors
    :param metric_min: float that sets the minimum value of the color scale; default is 0.0
    :param min_color: string of a color to be used to label the minimum value; must be a color accepted by matplotlib.colors; default is 'lightgray'
    :param under_color: string of a color to be used to label values below the minimum (used to denote null results); must be a color accepted by matplotlib.colors; default is 'white'
    :param colorbar_label: string to be used as the label for the colorbar; use mathtex writing for units if necessary; default is 'Metric'

    OUTPUTS:
    A new vis state file and colorbar figure that can be used to visualize the metric data on the structure.
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors
    from scipy import interpolate

    # ----------------------------------------
    # SETTING COLORS
    # ----------------------------------------
    # white used for residues with metric value below metric_min
    color_rgba = colors.to_rgba(under_color)
    color_defs = {}
    color_defs[33] = "color change rgb 33 %.3f %.3f %.3f\n"%(color_rgba[0],color_rgba[1],color_rgba[2])
    # min_color used for the minimum metric value
    min_color_rgba = np.array(colors.to_rgba(min_color))
    # max_color used for the maximum metric value
    max_color_rgba = np.array(colors.to_rgba(max_color))
    # color difference
    max_min_diff = max_color_rgba - min_color_rgba

    # ----------------------------------------
    # CREATING/FILLING COLOR DICTIONARIES
    # ----------------------------------------
    possible_colorids = list(range(34,1057))
    nColorids = len(possible_colorids)
    cmap_positions = np.linspace(0,1,nColorids)
    
    cdict = {'red':[], 'green':[], 'blue':[]}
    for colorid in possible_colorids:
        idx = colorid - possible_colorids[0]
        thiscolor = np.abs(min_color_rgba + idx * max_min_diff/(nColorids-1))
        ### VMD Stuff
        color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
        ### MATPLOTLIB Stuff
        cdict['red'].append((cmap_positions[idx],thiscolor[0],thiscolor[0]))
        cdict['green'].append((cmap_positions[idx],thiscolor[1],thiscolor[1]))
        cdict['blue'].append((cmap_positions[idx],thiscolor[2],thiscolor[2]))

    # ----------------------------------------
    # WRITING THE VIS STATE FILE
    # ----------------------------------------
    with open(vis_state_file_name,'w') as W:
        ### starting lines
        W.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')
    
        ### setting colorids
        W.write('# setting colorid rgb values\n')
        W.write(color_defs[33])
        for i in possible_colorids:
            W.write(color_defs[i])
    
        ### prepping the molecule and reps
        W.write('mol new ' + pdb_file_name + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
        W.write('mol delrep 0 top\n')
        W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
        W.write('mol color Beta\n')
        W.write('mol selection {all}\n')
        W.write('mol material AOEdgy\n')
        W.write('mol addrep top\n')
        W.write('mol scaleminmax top 0 %f %f\n'%(metric_min,metric_max))    # applies the correct colorbar scale to the newcartoon rep; must come after the "mol addrep top" line
        W.write('mol representation Licorice 0.150000 85.000000 85.000000\n')
        W.write('mol color Name\n')
        W.write('mol selection { protein }\n')
        W.write('mol material AOEdgy\n')
        W.write('mol addrep top\n\n')
    
        W.write('### setting viewpoints\n#REPLACE BELOW WITH NEW VIEWPOINTS LINE\nset viewpoints([molinfo top]) {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{0 0 0 0} {0 0 0 0} {0 0 0 0} {0 0 0 1}} {{0 0 0 0} {0 0 0 0} {0 0 0 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\nset topmol [molinfo top]\n\nforeach v $viewplist { \n  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\nforeach v $fixedlist {\n  molinfo $v set fixed 1\n}\nunset viewplist\nunset fixedlist\n')
    
    print('Finished creating vis-state file', vis_state_file_name)

    # ----------------------------------------
    # MATPLOTLIB Stuff
    # ----------------------------------------
    # creating a colorbar figure to accompany the VMD vis state colorbar;
    # depending on your renderer, the colors may not match with matplotlib's colorbar
    metric_range = metric_max - metric_min
    cmap = mpl.colors.LinearSegmentedColormap('my_cmap',cdict,nColorids)
    cmap.set_under(under_color)
    fig, ax = plt.subplots(figsize=(2,8))
    fig.subplots_adjust(right=0.5)
    norm = mpl.colors.Normalize(vmin=metric_min,vmax=metric_max)
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,extend='min',spacing='uniform',orientation='vertical',norm=norm,ticks=[metric_min,0.25*metric_range,0.50*metric_range,0.75*metric_range,metric_max])
    cb.set_label(r'%s'%(colorbar_label),size=16)
    cb.set_ticklabels([str(metric_min),str(0.25*metric_range),str(0.50*metric_range),str(0.75*metric_range),str(metric_max)])
    plt.savefig(colorbar_file_name,dpi=600,transparent=True)
    plt.close()
    
    print('Finished creating colorbar figure', colorbar_file_name)

    return


