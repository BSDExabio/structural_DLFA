
import sys
import argparse
import time
import platform
from pathlib import Path
import logging
import csv
import pandas
import pickle

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

import apoc_parsers
import rcsb_query
import uniprot_query

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


def append_timings(csv_writer, file_object, hostname, worker_id, start_time, stop_time,
                   query, return_code):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param file_object: file object associated with the CSV
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param query: query used as input to a task
    :param return_code: result of the task; 1=successful, 0=failure
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'query'      : query,
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


def _parse_alignment_score_file(alignment_results_file, file_type = 'TMalign'):
    """
    """
    start_time = time.time()
    worker = get_worker()
    aln_scores = []
    return_code = 0

    if file_type.upper() = 'APOC':
        aln_scores = apoc_parsers.parse_apoc_score_file(alignment_results_file)
        return_code = 1
    elif file_type.upper() = 'TMALIGN':
        pass
        #aln_scores = apoc_parsers.parse_apoc_score_file(alignment_results_file)

    return 1, aln_scores, platform.node(), worker.id, start_time, time.time(), return_code


def _query_rcsb(pdbid_chainid):
    """
    """
    start_time = time.time()
    worker = get_worker()
    uniprotid = ''
    return_code = 0
    try:
        uniprotid = rcsb_query.query_uniprot_str(pdbid_chainid)
        return_code = 1

    return 2, {pdbid_chainid: uniprotid}, platform.node(), worker.id, start_time, time.time(), return_code


def _query_uniprot_flat_file(uniprotid):
    """
    """
    start_time = time.time()
    worker = get_worker()
    meta_dict = {}
    return_code = 0
    try:
        meta_dict = uniprot_query.request_uniprot_metadata(uniprotid)
        return_code = 1

    return 3, {uniprotid: meta_dict}, platform.node(), worker.id, start_time, time.time(), return_code


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Parsing and requesting metadata associated with structural alignment hits.')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-list-file', '-inp', required=True, help='list file that contains the paths to alignment score files')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--tmscore-threshold', '-cut', required=True, help='float value between 0 and 1 used as the cutoff for TMscore results.')
    args = parser.parse_args()

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger',args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Alignment score files are listed in {args.input_list_file}')
    main_logger.info(f'TMscore value cutoff: {args.tmscore_threshold}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += f'Client information: {client}\n'
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')
    
    with open(args.input_list_file,'r') as alignment_list_file:
        alignment_files = alignment_list_file.readlines()
    main_logger.info(f'{len(alignment_files)} alignment results will be parsed.')

    # set up timing log file.
    timings_file = open(args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','query_pdb','task_type','return_code'])
    timings_csv.writeheader()

    # do the thing; step 1.
    step1_futures = client.map(_parse_alignment_score_file, alignment_files)

    # gather results from step1, step 2, and step3 into respective containers
    pdbid_chainid_list   = []
    uniprotid_list       = []
    protID_dict          = {}
    pdbid_to_uniprot_dict= {}
    uniprot_metadata_dict= {}
    ac = as_completed(step1_futures)
    for finished_task in ac:
        task_num, results, platform, workerid, start, stop, return_code = finished_task.result()
        # handling step 1 results:
        if task_num == 1:
            protID = alignment_pd_df['protein'][0]
            append_timings(timings_csv,timings_file,platform,workerid,start,stop,protID,return_code)
            main_logger.info(f'Structural alignment file associated with {protID} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
            protID_dict[protID] = alignment_pd_df.drop(labels['protein','Description'],axis=1) # removing redundant/unnecessary fields
            # loop over all pdbid_chainid's marked as "good" hits
            for pdbid_chainid in alignment_pd_df.log[alignment_pd_df['mscore'] > args.tmscore_threshold]['tname']:
                # only submit a new step 2 task if the pdbid_chainid has not been seen before
                if pdbid_chainid not in pdbid_chainid_list:
                    pdbid_chainid_list.append(pdbid_chainid)
                    ac.add(client.submit(_query_rcsb, pdbid_chainid))

        # handling step 2 results:
        if task_num == 2:
            pdbid_chainid = list(results.keys())[0]
            pdbid_to_uniprot_dict.update(results)
            append_timings(timings_csv,timings_file,platform,workerid,start,stop,pdbid_chainid,return_code)
            main_logger.info(f'The UniProt accession ID associated with {pdbid_chainid} has been queried. Return code: {return_code}. Took {stop-start} seconds.')
            # only submit a new step 3 task if the uniprot accession id != None or '' and also has not already been seen.
            if results[pdbid_chainid] and results[pdbid_chainid] not in uniprotid_list:
                uniprotid_list.append(results[pdbid_chainid])
                ac.add(client.submit(_query_uniprot_flat_file,results[pdbid_chainid])
            
        # handling step 3 results:
        if task_num == 3:
            uniprotid = list(results.keys())[0]
            uniprot_metadata_dict.update(results)
            append_timings(timings_csv,timings_file,platform,workerid,start,stop,uniprotid,return_code)
            main_logger.info(f'The flat file associated with {uniprotid} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
    
    # save dictionary of panda dataframes

    # save dictionary of the pdbid_chainid to uniprotid mapping

    # save dictionary of uniprot accession id meta data

    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

