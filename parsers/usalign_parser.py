#!/usr/bin/env python3
"""


"""
import time
from pathlib import Path
import pandas as pd

def parse_tmalign_score_file(fn):
    """UNNEEDED 
    """
    # TODO IMPLEMENT
    pass


def parse_tmalign_align_file(fn,alignment_type):
    """ Read given USalign alignment file, assuming alignment_type style/format

    :param fn: path string for score file
    :param alignment_type: string parameter to set format of alignment output to parse
    :return: pandas dataframe
    """
    start_time = time.time()
    results = {}
    with open('fn','r') as aln_log_file:
        lines = aln_log_files.readlines()

    if alignment_type.upper() == 'CP':
        # gather relevant lines
        model1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
        model2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
        aln_lines    = [line.strip() for line in lines if 'Aligned length=' in line]
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]

        # parse for relevant information
        results['model1']   = model1_lines[0].split()[3][:-2] # gather the path to the mobile structure
        results['model2']   = model2_lines[0].split()[3][:-2] # gather the path to the target structure
        results['tmscore1'] = float(re.search(r"(?<=TM-score= )(\d*\.\d*)",model1_lines[2]).group(1))
        results['tmscore2'] = float(re.search(r"(?<=TM-score= )(\d*\.\d*)",model2_lines[2]).group(1))
        results['Len1']     = int(re.search(r"(?<=L=)(\d*)",model1_lines[2]).group(1))
        results['Len2']     = int(re.search(r"(?<=L=)(\d*)",model2_lines[2]).group(1))
        results['d0_1']     = float(re.search(r"(?<=d0=)(\d*\.\d*)",model1_lines[2]).group(1))
        results['d0_2']     = float(re.search(r"(?<=d0=)(\d*\.\d*)",model2_lines[2]).group(1))
        # gather the alignment mapping
        results['map_2_to_1'] = {}
        # dictionary of tuples; keys are target residue index with values being (mobile residue index, target resname, mobile resname)
        for line in map_lines:
            temp = line.split()
            results['map_2_to_1'].extend({temp[7] : (temp[3],temp[5],temp[1])})
    
    elif alignment_type.upper() in ['SNS','FNS']:
        # gather relevant lines
        model1_lines = [line.strip() for line in lines if 'Structure_1:' in line]
        model2_lines = [line.strip() for line in lines if 'Structure_2:' in line]
        aln_lines    = [line.strip() for line in lines if 'Aligned length=' in line]
        map_lines    = [line.strip() for line in lines if ' CA ' == line[:4]]

        # parse for relevant information
        results['model1']   = model1_lines[0].split()[3][:-2] # gather the path to the mobile structure
        results['model2']   = model2_lines[0].split()[3][:-2] # gather the path to the target structure
        results['tmscore1'] = float(re.search(r"(?<=TM-score= )(\d*\.\d*)",model1_lines[2]).group(1))
        results['tmscore2'] = float(re.search(r"(?<=TM-score= )(\d*\.\d*)",model2_lines[2]).group(1))
        results['Len1']     = int(re.search(r"(?<=L=)(\d*)",model1_lines[2]).group(1))
        results['Len2']     = int(re.search(r"(?<=L=)(\d*)",model2_lines[2]).group(1))
        results['d0_1']     = float(re.search(r"(?<=d0=)(\d*\.\d*)",model1_lines[2]).group(1))
        results['d0_2']     = float(re.search(r"(?<=d0=)(\d*\.\d*)",model2_lines[2]).group(1))
        # gather the alignment mapping
        results['map_2_to_1'] = {}
        # dictionary of tuples; keys are target residue index with values being (mobile residue index, target resname, mobile resname, distance)
        for line in map_lines:
            temp = line.split()
            results['map_2_to_1'].extend({temp[7] : (temp[3],temp[5],temp[1],float(temp[8]))})

    else:
        print(f'{alignment_type} not expected by parser function. Check format type and retry.')
        return {}, start_time, stop_time, 1

    stop_time = time.time()
    return results, start_time, stop_time, 0


if __name__ == '__main__':
    # test harness for parsers

    df = parse_tmalign_score_file('/Users/may/Projects/data/Alphafold/TMalign-APoc/WP_010939197.1/WP_010939197.1_sco.dat')

    pass
