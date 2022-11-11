#!/usr/bin/env python3
""" Task manager for the dask pipeline for running TMalign on a set of alignments, saving log and alignment files to subdirectories. 
    USAGE: 
        python3 dask_taskmgr.py [-h] --scheduler-file SCHEDULER_FILE --alignment-list-file INPUT_FILE --script-path /path/to/dir/script.py --timings-file TIMINGS_FILE.csv --tskmgr-log-file TSKMGR.log
    INPUT: 
        -h, --help      show this help message and exit
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --alignment-list-file INPUT_FILE, -inp INPUT_FILE
                        list file that contains the paths to query and target pdb files as well as directory within which files will be written in
        --script-path /path/to/dir/script.py, -sp /path/to/dir/script.py
                        full path to the script to be run within the subprocess
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv 
                        CSV file for task processing timings
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the workflow will be written
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
                   return_code):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param return_code: return code from the subprocess call
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'return_code': return_code})
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
    """
    worker = get_worker()
    start_time = time.time()
    results_dict = {}
    query_protein  = proteins[0]
    target_proteins= proteins[1]
    output_directory = proteins[2]

    if type(target_proteins) == str:
        target_proteins = [target_proteins]

    for target in target_proteins:
        temp_path = Path(output_directory) / Path(target).stem
        temp_path.mkdir(mode=0o777,parents=True,exist_ok=True)
        
        try:
            completed_process = subprocess.run(f'bash {script} {query_protein} {target} {str(temp_path)}',shell=True,capture_output=True,check=True)
            return_code = completed_process.returncode
            if return_code != 0:
                print(query_protein, target, str(temp_path), return_code, file=sys.stderr, flush=True)

        except CalledProcessError as e:
            print(query_protein, target, e, file=sys.stderr, flush=True)
            return_code = 1
        
    stop_time = time.time()
    return platform.node(), worker.id, start_time, stop_time, return_code


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Molecular dynamics simulation task manager')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--alignment-list-file', '-inp', required=True, help='list file that contains the paths to query and target pdb files as well as directory within which files will be written in')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    args = parser.parse_args()

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Alignments are listed in {args.alignment_list_file}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')
    # number of workers
    main_logger.info(f'Client information: {client}')
    
    # set up timing log file.
    main_logger.info(f'Opening the timing file.')
    timings_file = open(args.timings_file, 'w')
    timings_csv  = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','return_code'])
    timings_csv.writeheader()

    # parse the query_pdb_list_file
    main_logger.info(f'Reading the alignment file.')
    with open(args.alignment_list_file,'r') as structures_file:
        alignments = [line.strip().split() for line in structures_file.readlines() if line[0] != '#']
    main_logger.info(f'Preparing to run {len(alignments):,} alignments.')
        
    # do the thing.
    aln_futures = client.map(submit_pipeline, alignments, script = args.script_path, pure=False)
    
    # gather alignment results.
    aln_completed = as_completed(aln_futures)
    for finished_task in aln_completed:
        hostname, worker_id, start_time, stop_time, return_code = finished_task.result()
        append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, return_code)
    
    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)
    #client.shutdown()

