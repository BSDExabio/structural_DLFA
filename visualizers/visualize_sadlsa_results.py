#!/usr/bin/env python3

"""
    Functions for pulling information out of SAdLSA pandas dataframes and visualizing results.
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


def ngl_colorscale(min_color, max_color, metric_max, metric_min, nColors = 1024, under_color = None, over_color = None):
    """
    :param     min_color:
    :param     max_color:
    :param   under_color:
    :param    over_color:
    :param       nColors: integer; default 1024 colors to be scaled through
    :return: number of colors in the scale, int
    :return: a list of hex codes
    :return: a matplotlib plot object
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors

    # ----------------------------------------
    # SETTING COLORS
    # ----------------------------------------
    color_defs = {}
    # set under and over colors
    if under_color:
        color_defs['under'] = colors.to_hex(under_color)
    if over_color:
        color_defs['over']  = colors.to_hex(over_color)
    
    # min_color used for the minimum metric value; need rgb color to perform the scaling
    min_color_rgba = np.array(colors.to_rgba(min_color))
    # max_color used for the maximum metric value; need rgb color to perform the scaling
    max_color_rgba = np.array(colors.to_rgba(max_color))
    # color difference
    max_min_diff = max_color_rgba - min_color_rgba

    # ----------------------------------------
    # CREATING/FILLING COLOR DICTIONARIES
    # ----------------------------------------
    cmap_positions = np.linspace(0,1,nColorids) # equally spaced colors within the range
    for colorid in range(nColors):
        thiscolor = np.abs(min_color_rgba + colorid * max_min_diff/(nColorids-1))   # calculate the rgb values for this iteration of color
        color_defs[colorid] = colors.to_hex(thiscolor)
        ### MATPLOTLIB Stuff
        cdict['red'].append((  cmap_positions[colorid],thiscolor[0],thiscolor[0]))
        cdict['green'].append((cmap_positions[colorid],thiscolor[1],thiscolor[1]))
        cdict['blue'].append(( cmap_positions[colorid],thiscolor[2],thiscolor[2]))
    
    # ----------------------------------------
    # MATPLOTLIB Stuff
    # ----------------------------------------
    # creating a colorbar figure to accompany the VMD vis state colorbar;
    # depending on your renderer, the colors may not match with matplotlib's colorbar
    metric_range = metric_max - metric_min
    cmap = colors.LinearSegmentedColormap('my_cmap',cdict,nColorids)
    if under_color:
        cmap.set_under(color_defs['under'])
    fig, ax = plt.subplots(figsize=(2,8))
    fig.subplots_adjust(right=0.5)
    norm = mpl.colors.Normalize(vmin=metric_min,vmax=metric_max)
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,extend='min',spacing='uniform',orientation='vertical',norm=norm,ticks=[metric_min,0.25*metric_range,0.50*metric_range,0.75*metric_range,metric_max])
    cb.set_label(r'%s'%(colorbar_label),size=16)
    cb.set_ticklabels([metric_min,'%.3f'%(0.25*metric_range),'%.3f'%(0.50*metric_range),'%.3f'%(0.75*metric_range),'%.3f'%(metric_max)])
    plt.savefig(colorbar_file_name,dpi=600,transparent=True)
    plt.close()
    
    print('Finished creating colorbar figure', colorbar_file_name)

    return nColors, color_defs, cdict


def vmd_colorscale(min_color, max_color, under_color = None, over_color = None):
    """
    :param     min_color:
    :param     max_color:
    :param   under_color:
    :param    over_color:
    :return: a list of colorids
    :return: a dictionary of r, g, and b values
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors

    # ----------------------------------------
    # SETTING COLORS
    # ----------------------------------------
    color_defs = {}
    min_color_idx = 33
    max_color_idx = 1057
    # set under and over colors
    if under_color:
        color_rgba = colors.to_rgba(under_color)
        # colorid 33 is the first color index of the color bar used in VMD
        color_defs[33] = "color change rgb 33 %.3f %.3f %.3f\n"%(color_rgba[0],color_rgba[1],color_rgba[2])
        min_color_idx = 34
    if over_color:
        color_rgba = colors.to_rgba(over_color)
        # colorid 1056 is the last color index of the color bar used in VMD
        color_defs[1056] = "color change rgb 33 %.3f %.3f %.3f\n"%(color_rgba[0],color_rgba[1],color_rgba[2])
        max_color_idx = 1056
    
    # min_color used for the minimum metric value
    min_color_rgba = np.array(colors.to_rgba(min_color))
    # max_color used for the maximum metric value
    max_color_rgba = np.array(colors.to_rgba(max_color))
    # color difference
    max_min_diff = max_color_rgba - min_color_rgba

    # ----------------------------------------
    # CREATING/FILLING COLOR DICTIONARIES
    # ----------------------------------------
    possible_colorids = list(range(min_color_idx,max_color_idx))
    nColorids = len(possible_colorids)
    cmap_positions = np.linspace(0,1,nColorids) # equally spaced colors within the range
    
    cdict = {'red':[], 'green':[], 'blue':[]}
    for colorid in possible_colorids:
        idx = colorid - possible_colorids[0]
        thiscolor = np.abs(min_color_rgba + idx * max_min_diff/(nColorids-1))
        ### VMD Stuff
        color_defs[colorid] = "color change rgb " + str(colorid) + " " + str(thiscolor[0]) + " " + str(thiscolor[1]) + " " + str(thiscolor[2]) + '\n'
        ### MATPLOTLIB Stuff
    if under_color:
        ngl_color_string += 'else if (atom.'+metric_string+' < '+str(metric_min)+') {return '+colors.to_hex(under_color)+'}'
    if over_color:
        ngl_color_string += 'else if (atom.'+metric_string+' > '+str(metric_min)+') {return '+colors.to_hex(over_color)+'}'
    ngl_color_string += end_string

    # ----------------------------------------
    # MATPLOTLIB Stuff
    # ----------------------------------------
    # creating a colorbar figure to accompany the VMD vis state colorbar;
    # depending on your renderer, the colors may not match with matplotlib's colorbar
    cmap = colors.LinearSegmentedColormap('my_cmap',cdict,nColors)
    metric_range = metric_max - metric_min

    if under_color:
        cmap.set_under(colors.to_hex(under_color))
    if over_color:
        cmap.set_over(colors.to_hex(over_color))

    fig, ax = plt.subplots(figsize=(2,8))
    fig.subplots_adjust(right=0.5)
    norm = colors.Normalize(vmin=metric_min,vmax=metric_max)


    if over_color and under_color:
        cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,extend='both',spacing='uniform',orientation='vertical',norm=norm)
    elif under_color:
        cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,extend='min',spacing='uniform',orientation='vertical',norm=norm)
    elif over_color:
        cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,extend='max',spacing='uniform',orientation='vertical',norm=norm)
    else:
        cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,spacing='uniform',orientation='vertical',norm=norm)


    cb.set_label(r'%s'%(colorbar_label),size=16)
    if nColors < 4:
        cb.set_ticks([metric_min,0.50*metric_range,metric_max])
        cb.set_ticklabels([metric_min,'%.3f'%(0.50*metric_range+metric_min),'%.3f'%(metric_max)])
    else:
        cb.set_ticks([metric_min,0.25*metric_range+metric_min,0.50*metric_range+metric_min,0.75*metric_range+metric_min,metric_max])
        cb.set_ticklabels([metric_min,'%.3f'%(0.25*metric_range),'%.3f'%(0.50*metric_range),'%.3f'%(0.75*metric_range),'%.3f'%(metric_max)])
    #plt.savefig(colorbar_file_name,dpi=600,transparent=True)
    #plt.close()
    plt.show()
    
    print('Finished creating colorbar figure')

    return nColors, ngl_color_string


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
    cb.set_ticklabels([metric_min,'%.3f'%(0.25*metric_range),'%.3f'%(0.50*metric_range),'%.3f'%(0.75*metric_range),'%.3f'%(metric_max)])
    plt.savefig(colorbar_file_name,dpi=600,transparent=True)
    plt.close()
    
    print('Finished creating colorbar figure', colorbar_file_name)

    return


