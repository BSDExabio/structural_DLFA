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


def grab_structure(pdbid,working_dir='./'):
    """
    Pull the structure from the RCSB database.
    
    INPUT:
    :param pdbid: string; 4 character long code associated with a structure on the RCSB
    :param working_dir: string; local or global path string pointing to where the files should be saved; option. TEMPORARY or NOTE: may not be maintained; hopefully we don't store any files long-term if at all.
    """

    from urllib import request

    pdb_urls = 'https://files.rcsb.org/view/%s.pdb' # static url for the pdb file associated with pdb ids; %s is replaced with the pdb id.
    fasta_urls = 'https://www.rcsb.org/pdb/download/viewFastaFiles.do?structureIdList=%s&compressionType=uncompressed' # static url for the fasta file associated with pdb ids; %s is replaced with the pdb id.
    cif_urls = 'https://files.rcsb.org/view/%s.cif' # static url for the mmcif file assocaited with pdb ids; %s is replaced with the pdb id.
    
    urls = [pdb_urls]
    #urls = [pdb_urls,fasta_urls,cif_urls]
    for url in urls:
        if url[-3:] == 'pdb':
            extension = 'pdb'
        elif url[-3:] == 'cif':
            extension = 'cif'
        elif "viewFastaFiles" in url:
            extension = 'fasta'
        else:
            print('something wrong is going on, %s'%(pdbid))
            return

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
        with open(working_dir + '/' + pdbid + '.' + extension,'w') as W:
            for line in response_list:
                W.write(line+'\n')

def edit_pdb(atom_data, metric_data, metric_type, string_type='file_path', default_value = -1.00):
    """
    Editing the atom data of a pdb-formatted data set; filling the 61-66 columns of the pdb to the float values of the metric_data.

    INPUTS:
    :param atom_data: kinda up in the air still. either a string of a global or local path to a PDB file OR the string object that contains all of the atomic data. 
    :param metric_data: a 2d array the contains the metric of interest; first column corresponds to the resids; 2nd column corresponds to the float values associated with the metric of interest. 
    :param metric_type: a string that acts as a descriptor for the metric; should not contain any spaces. TEMPORARY  
    
    OUTPUTS: 
    A new .pdb file that contains the metric of interest in the b-factor column. 
    """

    ### MAJOR MAPPING ERROR BETWEEN SEQUENCE POSITION REPORTED IN SADLSA ALIGNMENT RESULTS AND THE RESID COLUMN OF THE PDB... AAAAAAHHHHHHHHHH

    # starting to figure this out...

    import MDAnalysis
    file_name = atom_data + '_' + metric_type + '.pdb'
    u = MDAnalysis.Universe(atom_data)
    prot = u.select_atoms('protein and chainID %s'%(atom_data[-1])
    _all = u.select_atoms('chainID %s'%(atom_data[-1])
    _all.translate(-prot.center_of_mass())
    nRes = prot.n_residues
    for i in range(nRes):
        if prot.residues[i].resid in metric_data[:,0]:
            idx = np.argwhere(prot.residues[i].resid == metric_data[:,0])
            prot.residues[i].tempfactors = float(metric_data[idx,1])   # will be printed with 2 decimal points; 
        else:
            prot.residues[i].tempfactors = default_value
    prot.write(file_name)


def create_vmd_vis_state(vis_state_file_name, pdb_file_name, max_color, min_color='lightgray'):
    """
    """
    from scipy import interpolate
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors

    # ----------------------------------------
    # SETTING COLORS
    # ----------------------------------------
    # white used for residues with no metric value
    color_rgba = colors.to_rgba('white')
    color_defs[33] = "color change rgb 33 %.3f %.3f %.3f\n"%(color_rgba[0],color_rgba[1],color_rgba[2])
    # lightgray used for the minimum metric value
    min_color_rgba = np.array(colors.to_rgba('lightgray'))
    # user-set color for the maximum metric value
    max_color_rgba = np.array(colors.to_rgba(color))
    # color difference
    max_min_diff = max_color_rgba - min_color_rgba

    # ----------------------------------------
    # CREATING/FILLING COLOR DICTIONARIES
    # ----------------------------------------
    node_possible_colorids = list(range(34,1057))
    nNode_colorids = len(node_possible_colorids)
    cmap_positions = np.linspace(0,1,nNode_colorids)
    
    cdict = {'red':[], 'green':[], 'blue':[]}
    for colorid in node_possible_colorids:
        idx = colorid - node_possible_colorids[0]
        thiscolor = np.abs(min_color_rgba + idx * max_min_diff/(nNode_colorids-1))
        ### VMD Stuff
        color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
        ### MATPLOTLIB Stuff
        cdict['red'].append((cmap_positions[index],thiscolor[0],thiscolor[0]))
        cdict['green'].append((cmap_positions[index],thiscolor[1],thiscolor[1]))
        cdict['blue'].append((cmap_positions[index],thiscolor[2],thiscolor[2]))

    with open(vis_state_file_name,'w') as W:
        ### starting lines
        W.write('#!/usr/local/bin/vmd\nset viewplist {}\nset fixedlist {}\n\n')
    
        ### setting colorids
        W.write('# setting colorid rgb values\n')
        W.write(color_defs[33])
        for i in node_possible_colorids:
            W.write(color_defs[i])
    
        ### prepping the molecule and reps
        W.write('mol new ' + pdb_file_name + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
        W.write('mol delrep 0 top\n')
        W.write('mol representation NewCartoon 0.160000 50.000000 4.100000 0\n')
        W.write('mol color Beta\n')
        W.write('mol selection {all}\n')
        W.write('mol material AOEdgy\n')
        W.write('mol addrep top\n')
        W.write('mol representation Licorice 0.150000 85.000000 85.000000\n')
        W.write('mol color Name\n')
        W.write('mol selection { protein }\n')
        W.write('mol material AOEdgy\n')
        W.write('mol addrep top\n\n')
    
        W.write('### setting viewpoints\n#REPLACE BELOW WITH NEW VIEWPOINTS LINE\nset viewpoints([molinfo top]) {{{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}} {{0 0 0 0} {0 0 0 0} {0 0 0 0} {0 0 0 1}} {{0 0 0 0} {0 0 0 0} {0 0 0 0} {0 0 0 1}} {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}\nlappend viewplist [molinfo top]\nset topmol [molinfo top]\n\nforeach v $viewplist { \n  molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)\n}\nforeach v $fixedlist {\n  molinfo $v set fixed 1\n}\nunset viewplist\nunset fixedlist\n')
    
    print('Finished creating vis-state file', vis_state_file_name)
    return

