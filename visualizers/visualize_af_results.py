#!/usr/bin/env python3

"""
    For pulling information out of AF pandas dataframes. 
"""

def query_af_df(af_df,model_num):
    """
    Pull a specific column of data out of a pandas dataframe containing AlphaFold's pLDDT results for numerous models

    :param af_df: pandas dataframe object; filled with pLDDT results for numerous AF models
    :param model_num: string associated with columns in the dataframe; usd to pull a specific model's data
    :return: 1d pandas series object of a single model's pLDDT array; array indices are equivalent to the zero-indexed residue indices of the model. 
    """
    return af_df[model_num]


def create_vmd_vis_state(vis_state_file_name, colorbar_file_name, pdb_file_names_list, metric_max, max_color, metric_min=0., min_color='lightgray', under_color='white',colorbar_label='Metric'):
    """
    Writing the VMD visualization state file to view the aligned protein structures with the metric values in the b-factor columns of the pdbs
    Also, creating a colorbar figure that mirrors the color range used in the VMD vis state. 

    INPUTS:
    :param vis_state_file_name: string of a global or local path that the vis state file will be saved as
    :param colorbar_file_name: string of a global or local path that the colorbar will be saved as
    :param pdb_file_names_list: list of strings of a global or local path to the .pdb files that will be visualized
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
   
        # loop through pdbs, load them, and prep their representations. 
        for pdb in pdb_file_names_list:
            W.write('mol new ' + pdb + ' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n')
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


