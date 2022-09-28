
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
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, stop_time,
                   query, task_type, return_code):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param file_object: file object associated with the CSV
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param query: query used as input to a task
    :param task_type: integer used to denote which step of the workflow has been performed
    :param return_code: result of the task; 1=successful, 0=failure
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'query'      : query,
                         'task_type'  : task_type,
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


def _parse_alignment_score_file(alignment_results_file, file_type = 'APOC'):
    """
    """
    start_time = time.time()
    worker = get_worker()
    aln_scores_df = []
    return_code = 0

    if file_type.upper() == 'APOC':
        aln_scores_df = apoc_parsers.parse_apoc_score_file(alignment_results_file)
        return_code = 1
    elif file_type.upper() == 'TMALIGN':
        pass
        #aln_scores_df = apoc_parsers.parse_apoc_score_file(alignment_results_file)

    return 1, aln_scores_df, platform.node(), worker.id, start_time, time.time(), return_code


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
    except:
        print(f'failed to pull the uniprot accession id associated with {pdbid_chainid}. oh well...')

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
    except:
        print(f'failed to pull {uniprotid} flat file. oh well...')

    return 3, {uniprotid: meta_dict}, platform.node(), worker.id, start_time, time.time(), return_code


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Parsing and requesting metadata associated with structural alignment hits.')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-list-file', '-inp', required=True, help='list file that contains the paths to alignment score files. Stem of this file name will be used to name the parsed alignment results file.')
    parser.add_argument('--output-dir', '-out', required=True, help='path to a directory (already made) where all output files will be stored.')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    parser.add_argument('--tmscore-threshold', '-cut', required=True, help='float value between 0 and 1 used as the cutoff for TMscore results.')
    #parser.add_argument('--uniprot-metadata-pickle', '-meta', help="path to a pickle file associated with the uniprot accession IDs' metadata dictionary.")
    #parser.add_argument('--pdbid-to-uniprot-map-pickle', '-map', help='path to a pickle file associated with the mapping between pdbids and uniprot accession ids.')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    tmscore_threshold = float(args.tmscore_threshold)
    
    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger', str(output_dir / args.tskmgr_log_file))
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    main_logger.info(f'Alignment score files are listed in {args.input_list_file}')
    main_logger.info(f'TMscore value cutoff: {tmscore_threshold}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += f'Client information: {client}\n'
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')
    
    with open(args.input_list_file,'r') as alignment_list_file:
        alignment_files = [line.strip() for line in alignment_list_file.readlines()]
    main_logger.info(f'{len(alignment_files)} alignment results will be parsed.')
    stem = Path(args.input_list_file).stem

    # set up timing log file.
    timings_file = open( str(output_dir / args.timings_file), 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','query','task_type','return_code'])
    timings_csv.writeheader()

    ### do the thing; only submitting step 1 for now since we know we need to run this set of tasks.
    step1_futures = client.map(_parse_alignment_score_file, alignment_files)

    ### gather results from step1, step 2, and step3 into respective containers
    protID_dict          = {}   # home for this run's parsed structural alignment results
    
    #NOTE... add code here
    # if provided as an argument, read in a map file of pdbid_to_uniprot; 
    # update this map dictionary object as new pdbid_to_uniprotid runs are performed
    # avoid redundant collection of metadata associated with a pdbid if we've already seen't it
    pdbid_chainid_list   = []
    pdbid_to_uniprot_dict= {}
    
    #NOTE... add code here
    # if provided as an argument, read in a dictionary filled with uniprot accession id metadata 
    # update this dictionary object as new uniprot accession ids are seen and parsed.
    uniprotid_list       = []
    uniprot_metadata_dict= {}

    # collect futures into a bucket that we will add tasks to
    ac = as_completed(step1_futures)
    for finished_task in ac:
        task_num, results, hostname, workerid, start, stop, return_code = finished_task.result()
        # handling step 1 results:
        if task_num == 1:
            protID = results['protein'][0]
            append_timings(timings_csv,timings_file,hostname,workerid,start,stop,protID,task_num,return_code)
            main_logger.info(f'Structural alignment file associated with {protID} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
            if return_code == 1:
                # removing redundant/unnecessary fields from the pandas dataframe object; store away in the protID_dict dictionary object
                protID_dict[protID] = results.drop(labels=['protein','Description'],axis=1)
                # loop over all pdbid_chainid's marked as "good" hits ('mscore' > tmscore_threshold)
                for pdbid_chainid in results.loc[results['mscore'] > tmscore_threshold]['tname']:
                    # only submit a step 2 task if the pdbid_chainid has never been seen before
                    if pdbid_chainid not in pdbid_chainid_list:
                        pdbid_chainid_list.append(pdbid_chainid)
                        ac.add(client.submit(_query_rcsb, pdbid_chainid))

        # handling step 2 results:
        elif task_num == 2:
            pdbid_chainid = list(results.keys())[0]
            pdbid_to_uniprot_dict.update(results)
            append_timings(timings_csv,timings_file,hostname,workerid,start,stop,pdbid_chainid,task_num,return_code)
            main_logger.info(f'The UniProt accession ID associated with {pdbid_chainid} has been queried. Return code: {return_code}. Took {stop-start} seconds.')
            # only submit a new step 3 task if the return_code is 1, the uniprot accession id != None or '', and the uniprot aaccession id has not already been seen.
            if return_code == 1 and results[pdbid_chainid] and results[pdbid_chainid] not in uniprotid_list:
                uniprotid_list.append(results[pdbid_chainid])
                ac.add(client.submit(_query_uniprot_flat_file,results[pdbid_chainid]))
            
        # handling step 3 results:
        elif task_num == 3:
            uniprotid = list(results.keys())[0]
            uniprot_metadata_dict.update(results)
            append_timings(timings_csv,timings_file,hostname,workerid,start,stop,uniprotid,task_num,return_code)
            main_logger.info(f'The flat file associated with {uniprotid} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
    
    # save dictionary of panda dataframes
    with open( str(output_dir / f'{stem}_structural_alignment_results.pkl'), 'w') as out:
        pickle.dump(protID_dict,out,protocol=pickle.HIGHEST_PROTOCOL)

    # save dictionary of the pdbid_chainid to uniprotid mapping
    with open( str(output_dir / 'pdbid_to_uniprotid_map.pkl'), 'w') as out:
        pickle.dump(pdbid_to_uniprot_dict,out,protocol=pickle.HIGHEST_PROTOCOL)

    # save dictionary of uniprot accession id meta data
    with open( str(output_dir / 'uniprot_metadata.pkl'), 'w') as out:
        pickle.dump(uniprot_metadata_dict,out,protocol=pickle.HIGHEST_PROTOCOL)

    # gotta finish the annotation pipeline work... matching structure alignment hits to the uniprot metadata... maybe for another script?

    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

