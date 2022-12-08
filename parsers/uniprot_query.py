#!/usr/bin/env python3
"""
    Request Uniprot flat files associated with given UniProt Accession IDs

    (Written by Russell B Davidson, davidsonrb@ornl.gov, rbdavid)
"""

import requests
import datetime
import copy

flat_file_url = 'https://www.uniprot.org/uniprot/%s.txt'
dictionary_namespace = {'entry_name': '',           # ID line
                        'status': '',               # ID line
                        'sequence_length': '',      # ID line
                        'original_date': '',        # 1st DT line
                        'sequence_date': '',        # 2nd DT line
                        'modified_date': '',        # 3rd DT line
                        'gene_name': '',            # GN line
                        'species_name': '',         # OS line
                        'database_references': [],  # DR lines
                        'primary_evidence': '',     # PE line
                        'features': [],             # FT lines
                        'sequence': '>',            # SQ and blank lines afterwards
                        'ecIDs': [],                # Gather any lines that have EC IDs in them
                        'status': 0}

def request_uniprot_metadata(accession_id):
    """
    """
    try:
        response = requests.get(flat_file_url %(accession_id))
    except Exception as e:
        print(f'Requesting the associated flat file for {accession_id} failed.')
        return {'status':-9999}
    
    status_code   = response.status_code
    response_text = response.text
   
    # successful response
    if status_code == 200:
        uni_dict = copy.deepcopy(dictionary_namespace)
        uni_dict['status'] = 0
        response_lines = response_text.split('\n')
        for line in response_lines:
            # gathering any and all instances where  EC ids are presented... morbid curiosity
            if 'EC=' in line:
                uni_dict['ecIDs'].append(line)
            # end of file is denoted by '//' so break from the for loop if it occurs
            if line == '//':    # used to denote end of file for the UniProt flat files
                break
            # parse first line
            elif line[:2] == 'ID':
                 temp = line.split()
                 uni_dict['entry_name'] = temp[1]
                 uni_dict['status'] = temp[2][:-1] # stupid semicolon at the end of the status string
                 uni_dict['sequence_length'] = int(temp[3])
            # parse the chronology lines
            elif line[:2] == 'DT':
                 if 'integrated' in line:
                     uni_dict['original_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
                 elif 'sequence version' in line:
                     uni_dict['sequence_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
                 elif 'entry version' in line:
                     uni_dict['modified_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
            # parse the gene name line(s)
            elif line[:2] == 'GN':
                uni_dict['gene_name'] += line[5:]
            # parse the species name line(s)
            elif line[:2] == 'OS':
                uni_dict['species_name'] += line[5:]
            # parse the database reference lines
            elif line[:2] == 'DR':
                uni_dict['database_references'].append([line[5:]])
            # parse the primary evidence line(s)
            elif line[:2] == 'PE':
                uni_dict['primary_evidence'] += line[5:]
            # parse the features lines
            # if the line is associated with a new feature, the 5th character
            # will be the first letter of the feature type and not a ' '
            # in this case, append a new list to the features list.
            elif line[:2] == 'FT' and line[5] != ' ':
                uni_dict['features'].append(line[5:].split())
            # if the line is associated with a previous feature, the 5th
            # character will be a ' '. 
            # in this case, append the line's information to the previous list
            elif line[:2] == 'FT' and line[5] == ' ':
                uni_dict['features'][-1].append(line[2:].strip())
            # parse the sequence lines
            # the 'SQ' key is the notification that the sequence information
            # begins on this line
            elif line[:2] == 'SQ' and line[5] != ' ':
                uni_dict['sequence'] += line[5:] + '\n'
            # the only instance where the first two characters of a line == '  '
            # is a line that continues with sequence information
            # of course, these lines are formatted; blocks of 10 AAs separated
            # by a space character; 60 AAs on a line
            elif line[:2] == '  ':
                temp = line.split()
                seq_string = ''.join(temp)
                uni_dict['sequence'] += seq_string
        
        return uni_dict

    # if the request response fails with known 'failure' status codes
    elif status_code in [400,404,410,500,503]:
        print(uniprotID, status_code, response_text)
        return {'status':status_code}
    # if the request response fails for any other reason
    else:
        print(uniprotID, status_code, 'Something really funky is going on')
        return {'status':status_code}


if __name__ == '__main__':
    # real uniprot accession id
    results = request_uniprot_metadata('A0A0M3KL33')
    assert results['status'] == 0

    # fake uniprot accession id
    results = request_uniprot_metadata('BOGUS')
    assert results['status'] != 0

