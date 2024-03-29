#!/usr/bin/env python3
""" Task manager for running the dask pipeline for running TMalign of a set of query protein structures against a library of pdb structures.
    USAGE: 
        python3 dask_taskmgr.py [-h] --scheduler-file SCHEDULER_FILE --query-pdb-list-file INPUT_FILE --target-pdb-list-file TARGET_FILE --timings-file TIMINGS_FILE.csv --working-dir /path/to/dir/ --script-path /path/to/dir/script.py --tskmgr-log-file TSKMGR.log
    INPUT: 
        -h, --help      show this help message and exit
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --query-pdb-list-file INPUT_FILE, -inp INPUT_FILE
                        list file that contains the paths and sequence lengths of protein models
        --target-pdb-list-file LIB_DIR_PATH, -lib LIB_DIR_PATH
                        path to a directory within which a library of protein models are held
        --script-path /path/to/dir/script.py, -sp /path/to/dir/script.py
                        full path to the script to be run within the subprocess
        --working-dir /path/to/dir/, -wd /path/to/dir/
                        full path to the directory within which files will be written
        --subdirectory-string DIRNAME, -sub DIRNAME
                        string to be used for consistent naming of subdirectories within which results files will be saved for each query structure.
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv 
                        CSV file for task processing timings
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the workflow will be written
        --nRanked-structures INTEGER_VALUE, -nrank INTEGER_VALUE
                        integer value set for number of top scoring alignments hits to be written to ranked file
        --sorted-by COLUMN_NAME, -sort COLUMN_NAME
                        string to set which column of the pandas dataframe to rank the alignment hits with. options: 'maxTMscore', 'TMscore1', 'TMscore2', ...
        --num-workers nWORKERS_IN_WORKFLOW, -nw nWORKERS_IN_WORKFLOW
                        integer value used to set the number of workers available to be used
"""

import time
import argparse
import platform
import numpy as np
from pathlib import Path
import sys
import logging
import itertools
import pickle
import pandas
import csv

import subprocess
from subprocess import CalledProcessError

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

#######################################
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, stop_time,
                   nSuccesses):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param nSuccesses: result of the subprocess call
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'nSuccesses' : nSuccesses})
    file_object.flush()


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


#######################################
### DASK RELATED FUNCTIONS
#######################################

def get_num_workers(client):
    """ Get the number of active workers
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


def submit_pipeline(proteins, script):
    """
    function used in the dask workflow to define a task
    performs the alignments between proteins[0] (a string) and proteins[1] (a list)
    the type of alignment performed is controlled by the script parameter; this is interchangeable as long as output from the script follows the _expected format_ (see below)
    
    EXPECTED FORMAT:
    ```
    #PDBchain1	PDBchain2	TM1	TM2	RMSD	ID1	ID2	IDali	L1	L2	Lali
    mobile.pdb	target.pdb	0.5970	0.5970	4.34	0.540	0.540	0.786	448	448	308
    Total CPU time is  0.64 seconds
    '''
    The 2nd line holds the relevant quantitative metrics we wish to gather.

    If format changes, this chunk of code will need to be updated to follow. 
    
    INPUT:
        :param proteins: list with elements 0 and 1 being a string and a list, respectively
                         proteins[0] is a path string that points to a pdb file used as the mobile structure
                         proteins[1] is a list of path strings that point to the pbd files to be used as the target structures
        :param script: path string that points to the bash script that will be subprocess.run'd to automate the alignment between mobile and target structures
    OUTPUT:
        :return: resource/node information, dask worker id, start time, stop time, number of alignments performed successfully, dictionary that houses the results
    """
    worker = get_worker()
    start_time = time.time()
    results_dict = {}
    query_protein  = proteins[0]
    target_proteins= proteins[1]
    if type(target_proteins) == str:
        target_proteins = [target_proteins]

    for target in target_proteins:
        try:
            completed_process = subprocess.run(f'bash {script} {query_protein} {target}',shell=True,capture_output=True,check=True)
            results_dict[f'{query_protein}_{target}'] = [float(elem) for elem in completed_process.stdout.decode().split('\n')[1].split('\t')[2:]]   # [TMscore1,TMscore2,RMSD,SeqID1,SeqID2,SeqIDali,L1,L2,Lali]
            #results_dict[f'{query_protein}_{target}'] = [float(elem) for elem in completed_process.stdout.split()]

        except CalledProcessError as e:
            print(query_protein, target, e, file=sys.stderr, flush=True)
        
    stop_time = time.time()
    return platform.node(), worker.id, start_time, stop_time, len(list(results_dict.keys())), results_dict


def post_analysis_pipeline(query_str, result_files_list, outputdir_str= './', subdir_str='TMalign', nRanked=100, sorted_by='maxTMscore'):
    """
    """
    start_time = time.time()
    protein_id = Path(query_str).parent.name
    temp_path = Path(outputdir_str) / protein_id / subdir_str
    temp_path.mkdir(mode=0o777,parents=True,exist_ok=True)
    # open and write alignment results to file
    with open(str(temp_path) + '/alignment_results.dat','w') as out_file, open(str(temp_path) + '/alignment_results.log','w') as log_file:
        out_file.write(f'Target Path,TMscore1,TMscore2,RMSD,SeqID1,SeqID2,SeqIDali,Len1,Len2,LenAligned\n') # [TMscore1,TMscore2,RMSD,SeqID1,SeqID2,SeqIDali,L1,L2,Lali]
        for result_file in result_files_list:
            with open(result_file,'rb') as infile:
                temp_results = pickle.load(infile)
            for key, value in temp_results.items():
                if query_str not in key:
                    continue
                elif len(value) != 9:
                    log_file.write( f"{key}, {value}, len(value) does not match expected length. something wrong with the TMalign results?\n")
                    continue
                else:
                    target = '/' + key.split('_/')[1]
                    out_file.write(f'{target},{value[0]},{value[1]},{value[2]},{value[3]},{value[4]},{value[5]},{value[6]},{value[7]},{value[8]}\n')
        
    # read in the results file and parse
    df = pandas.read_csv(str(temp_path) + '/alignment_results.dat')
    # create a maxTMscore column if expecting to sort by this value
    if sorted_by.upper() == 'MAXTMSCORE':
        df['maxTMscore'] = np.max(df[['TMscore1','TMscore2']],axis=1)
    # check that sorted_by is one of the column names of the panda dataframe
    elif sorted_by not in list(df.columns):
        print(f'{sort_by} not in pandas dataframe columns {list(df.columns)}. No ranked_alignment_results.dat file is written.', file=sys.stderr, flush=True)
        stop_time = time.time()
        return start_time, stop_time, query_str

    # sort the dataframe
    sorted_df = df.sort_values(sorted_by,ascending=False)
    # write the sorted dataframe to file
    try:
        sorted_df.head(nRanked).to_csv(str(temp_path) + '/ranked_alignment_results.dat',index=False)
    # if nRanked > the total number of alignments performed, the above code will raise an exception
    # instead, just write the full set of ranked alignment results. 
    except:
        sorted_df.to_csv(str(temp_path) + '/ranked_alignment_results.dat',index=False)

    stop_time = time.time()
    return start_time, stop_time, query_str


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Molecular dynamics simulation task manager')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--query-pdb-list-file', '-inp', required=True, help='list file that contains the paths and sequence lengths of protein models')
    parser.add_argument('--target-pdb-list-file', '-lib', required=True, help='list file that contains the paths and sequence lengths of protein models assocaited with the target library')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    parser.add_argument('--working-dir', '-wd', required=True, help='path that points to the working directory for the output files')
    parser.add_argument('--subdirectory-string', '-sub', required=True, help='string to be used for consistent naming of subdirectories within which results files will be saved for each query structure.')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--nRanked-structures', '-nrank', required=True, help='integer value set for number of top scoring alignment hits to ranked file')
    parser.add_argument('--sorted-by', '-sort', required=True, help='string used to denote which column to rank structures by')
    parser.add_argument('--num-workers', '-nw', required=True, help='integer used to set how many workers are being claimed by the client')
    args = parser.parse_args()

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Working directory: {args.working_dir}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Query models are listed in {args.query_pdb_list_file}')
    main_logger.info(f'Target models are listed in {args.target_pdb_list_file}')
    main_logger.info(f'Alignment results will be saved within subdirectories named {args.subdirectory_string}')
    main_logger.info(f'All alignment results will be ranked by the {args.sorted_by} metric')
    main_logger.info(f'{args.nRanked_structures} ranked alignment results will be saved')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')
    # number of workers
    NUM_WORKERS = int(args.num_workers)
    main_logger.info(f'Client information: {client} {NUM_WORKERS}')
    
    # set up timing log file.
    main_logger.info(f'Opening the timing file.')
    timings_file = open(args.timings_file, 'w')
    timings_csv  = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','nSuccesses'])
    timings_csv.writeheader()

    # parse the query_pdb_list_file
    main_logger.info(f'Reading the query structure file.')
    with open(args.query_pdb_list_file,'r') as structures_file:
        query_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    # sorted largest to smallest of the query/mobile structure list; system size is a basic but good estimate of computational cost
    sorted_query_list = sorted(query_list, key = lambda x: int(x[1]))[::-1]
    del query_list
    sorted_query_list = [elem[0] for elem in sorted_query_list]
    main_logger.info(f'Preparing to run {len(sorted_query_list):,} query structures, sorted by longest to shortest sequences.')

    # parse the target_pdb_list_file
    main_logger.info(f'Reading the target structure file.')
    with open(args.target_pdb_list_file,'r') as structures_file:
        target_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    
    # sorted largest to smallest of the target structure list; system size is a basic but good estimate of computational cost
    sorted_target_list = sorted(target_list, key = lambda x: int(x[1]))[::-1]
    del target_list
    sorted_target_list = [elem[0] for elem in sorted_target_list]
    nTargets = len(sorted_target_list)
    main_logger.info(f'Preparing to run alignments against {nTargets:,} target structures.')
    
    # determining total number of alignments to be perforemd
    nAlignments = len(sorted_query_list)*len(sorted_target_list)
    
    ############
    # if total number of alignment tasks is relatively small, 
    # then create the full iterable list and map directly to client
    ############
    if nAlignments < 10**4 or nTargets < NUM_WORKERS:
        main_logger.info(f'The total number of alignments ({nAlignments:,}) is relatively small. Running all alignments as individual tasks.')
        
        # do the thing.
        aln_futures = client.map(submit_pipeline, list(itertools.product(sorted_query_list,sorted_target_list)), script = args.script_path, pure=False) #, batch_size=10**4

        # gather alignment results.
        results_dict = {}
        aln_completed = as_completed(aln_futures)
        for finished_task in aln_completed:
            hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
            query_target = list(results.keys())[0]
            main_logger.info(f'{query_target} alignment has finished.')
            append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, query_target)
            results_dict.update(results)
        
        out = f'{args.working_dir}/results_dictionary.pkl'
        main_logger.info(f'RESULTS WRITTEN TO {out}. {time.time()}')
        with open(out,'wb') as out_file:
            pickle.dump(results_dict,out_file)

        # submitting the post_analysis_pipeline as a dask task for all query structures
        parse_futures = client.map(post_analysis_pipeline, sorted_query_list, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, nRanked = int(args.nRanked_structures), sorted_by = args.sorted_by, pure=False)
        parse_completed = as_completed(parse_futures)
        for finished_task in parse_completed:
            start_time, stop_time, query = finished_task.result()
            main_logger.info(f'Parsed alignment hits for {query}, taking {stop_time - start_time} seconds.')

    ############
    # if total number of alignments is not relatively small (>10k), 
    # then partition the target list and map query-target partitions to the client
    ############
    else:
        main_logger.info(f'The total number of alignments ({nAlignments:,}) is large. Will run the alignments in batches of targets.')
        target_sublists = [sorted_target_list[i::NUM_WORKERS*3] for i in range(NUM_WORKERS*3)]
        nTargets_per_sublist = int(np.mean([len(sublist) for sublist in target_sublists]))
        nTasks = len(sorted_query_list)*len(target_sublists)
        main_logger.info(f'The total number of tasks is {nTasks:,} where each task is the alignment of a query structure to a set of target structures (average lengths of {nTargets_per_sublist}).')
        for query in sorted_query_list:
            query_start = time.time()
            # do the thing.
            aln_futures = client.map(submit_pipeline, list(itertools.product([query],target_sublists)), script = args.script_path, pure=False) #, batch_size=10**4
            
            # gather query's alignment results.
            results_dict = {}
            aln_completed = as_completed(aln_futures)
            for finished_task in aln_completed:
                hostname, worker_id, start_time, stop_time, nSuccesses, results = finished_task.result()
                append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, nSuccesses)
                results_dict.update(results)
        
            out = f'{args.working_dir}/results.pkl'
            with open(out,'wb') as out_file:
                pickle.dump(results_dict,out_file)
    
            post_analysis_pipeline(query, result_files_list = [out], outputdir_str = args.working_dir, subdir_str = args.subdirectory_string, nRanked = int(args.nRanked_structures), sorted_by = args.sorted_by)
            main_logger.info(f'Finished running, collecting, and outputting alignment results for {query}. This took {time.time() - query_start} seconds.')
    
    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)
    #client.shutdown()
