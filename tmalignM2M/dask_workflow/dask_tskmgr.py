#!/usr/bin/env python3
""" Task manager for running the dask pipeline for running TMalign of a set of query protein structures against a library of pdb structures.

    USAGE: 
        python3 dask_taskmgr.py [-h] --scheduler-file SCHEDULER_FILE --input-pdb-list INPUT_FILE --library-pdb_directory LIB_DIR_PATH --timings-file TIMINGS_FILE.csv --working-dir /path/to/dir/ --script-path /path/to/dir/script.py --tskmgr-log-file TSKMGR.log

    INPUT: 
        -h, --help      show this help message and exit
        --scheduler-file SCHEDULER_FILE, -s SCHEDULER_FILE
                        dask scheduler file
        --input-pdb-list INPUT_FILE, -inp INPUT_FILE
                        list file that contains the paths and sequence lengths of protein models
        --tmscore-threshold THRESHOLD, -cut THRESHOLD
                        float value between 0 and 1 used as the cutoff for TMscore results.
        --library-pdb-directory LIB_DIR_PATH, -lib LIB_DIR_PATH
                        path to a directory within which a library of protein models are held
        --timings-file TIMINGS_FILE.csv, -ts TIMINGS_FILE.csv 
                        CSV file for task processing timings
        --working-dir /path/to/dir/, -wd /path/to/dir/
                        full path to the directory within which files will be written
        --script-path /path/to/dir/script.py, -sp /path/to/dir/script.py
                        full path to the script to be run within the subprocess
        --tskmgr-log-file TSKMGR.log, -log TSKMGR.log 
                        path for a file within which logging output for the workflow will be written

"""

import time
import argparse
import platform
import os
from pathlib import Path
import logging

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
                   return_code, query_pdb):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param return_code: result of the subprocess call
    :param query_pdb: model file that was processed
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'return_code': return_code,
                         'query_pdb'  : query_pdb})
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


def submit_pipeline(query_pdb, library_pdb_directory, script, working_directory, threshold):
    """
    """
    worker = get_worker()
    start_time = time.time()

    save_directory = Path(str(working_directory) + '/' + Path(query_pdb).parent.name + '/TMAlign')
    save_directory.mkdir(mode=0o777,parents=True,exist_ok=True)
    
    return_code = 1
    try:
        completed_process = subprocess.run(f'bash {script} {query_pdb} {library_pdb_directory} {str(save_directory)} {threshold}',shell=True,capture_output=True,check=True,cwd=working_directory)
        return_code = 0

    except CalledProcessError as e:
        print(str(e), file=sys.stderr, flush=True)

    finally:
        print(platform.node(), worker.id, start_time, time.time(), return_code, query_pdb)
        return platform.node(), worker.id, start_time, time.time(), return_code, query_pdb


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Molecular dynamics simulation task manager')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-pdb-list', '-inp', required=True, help='list file that contains the paths and sequence lengths of protein models')
    parser.add_argument('--library-pdb-directory', '-lib', required=True, help='path to a directory within which a library of protein models are held')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--working-dir', '-wd', required=True, help='path that points to the working directory for the output files')
    parser.add_argument('--script-path', '-sp', required=True, help='path that points to the script for the subprocess call')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--tmscore-threshold', '-cut', required=True, help='float value between 0 and 1 used as the cutoff for TMscore results.')
    args = parser.parse_args()

    if args.working_dir[-1] != os.path.sep:
        args.working_dir += os.path.sep

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Working directory: {args.working_dir}')
    main_logger.info(f'Path to subprocess script: {args.script_path}')
    main_logger.info(f'Query models are listed in {args.input_pdb_list}')
    main_logger.info(f'Library models are found in {args.library_pdb_directory}')
    main_logger.info(f'TMscore value cutoff for saving results: {args.tmscore_threshold}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')
    main_logger.info(f'Client information: {client}')
    NUM_WORKERS = get_num_workers(client)

    # parse the input_pdb_list
    with open(args.input_pdb_list,'r') as structures_file:
        structures_list = [line.split() for line in structures_file.readlines() if line[0] != '#']
    main_logger.info(f'{structures_list}')
    # sorted largest to smallest of the structure list; system size is a basic but good estimate of computational cost
    sorted_structures_list = sorted(structures_list, key = lambda x: int(x[1]))[::-1]
    main_logger.info(f'{sorted_structures_list}')
    sorted_structures_list = [elem[0] for elem in sorted_structures_list]
    main_logger.info(f'{sorted_structures_list}')
    main_logger.info(f'Preparing to run {len(sorted_structures_list)} tasks, sorted by longest to shortest sequences.')

    # set up timing log file.
    timings_file = open(args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','return_code','query_pdb'])
    timings_csv.writeheader()

    # do the thing.
    task_futures = client.map(submit_pipeline, sorted_structures_list, library_pdb_directory = args.library_pdb_directory, script = args.script_path, working_directory = args.working_dir, threshold = float(args.tmscore_threshold), pure=False)

    # gather results.
    ac = as_completed(task_futures)
    for finished_task in ac:
        hostname, worker_id, start_time, stop_time, return_code, query_pdb = finished_task.result()
        main_logger.info(f'{query_pdb} has finished, return code: {return_code}.')
        append_timings(timings_csv, timings_file, hostname, worker_id, start_time, stop_time, return_code, query_pdb)

    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)
    #client.shutdown()

