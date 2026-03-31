#!/usr/bin/env python3
import argparse
import os
import pathlib
import time
from collections import defaultdict
from multiprocessing import Manager, Process, Queue, Value

import arrow
import numpy as np
import pandas as pd
import pysam
from natsort import natsorted

import warnings
warnings.filterwarnings("ignore", category=FutureWarning, message=".*swapaxes.*")


def build_sites_index(sites_set):
    sites_index = {}
    for chromosome, chrom_sites in sites_set.groupby("chromosome", sort=False):
        sites_index[chromosome] = {
            "forward": chrom_sites.set_index("start")[["epic_id"]],
            "reverse": chrom_sites.set_index("end")[["epic_id"]],
        }
    return sites_index


def iter_unique_methylation_alignments(alignments_by_read):
    for alignments in alignments_by_read.values():
        if len(alignments) != 1:
            continue
        read_alignment = alignments[0]
        if read_alignment[2] == {}:
            continue
        yield read_alignment


def chunk_records(records, chunk_size):
    if chunk_size <= 0:
        chunk_size = 1
    for idx in range(0, len(records), chunk_size):
        yield records[idx : idx + chunk_size]


def read_bam_files(path, recursive, file_queue, bams_amount, bam_filter):
    # function to read BAM files from minimap2 output  

    # params: 
    # path: directory where BAM files are looked for 
    # recursive: boolean whether searching for files should be recursive
    # file_queue: queue that stores BAM files for further processing
    # bams_amount: counter for tracking total number of found BAM files 
    # bam_filter: optional string for filtering on barcodes 

    d = pathlib.Path(path)
    if recursive:
        # if no filter provided, queue all BAM files 
        if bam_filter is None:
            with bams_amount.get_lock():
                # use natsorted for natural file sorting and rglob for recursive globbing
                bams_amount.value = len(natsorted(d.rglob('*.bam')))
            # queue each file 
            for file in natsorted(d.rglob('*.bam')):
                file_queue.put(str(file))
        # if filter is provided, queue all files that contain the filter string 
        else:
            for file in natsorted(d.rglob('*.bam')):
                if bam_filter in str(file):
                    file_queue.put(str(file))
                    with bams_amount.get_lock():
                        bams_amount.value += 1
    else:
        with bams_amount.get_lock():
            bams_amount.value = len(natsorted(d.glob('*.bam')))
        # glob instead of rglob for non-recursive search     
        for file in natsorted(d.glob('*.bam')):
            file_queue.put(str(file))

def io_handler(
    file_queue,
    methylation_queue,
    methylation_read_number,
    bams_analysed,
    runs,
    read_amount,
    chunk_size,
):
    # function for processing of sequencing data from BAM files 

    # params: 
    # file_queue: queue of BAM files to process
    # methylation_queue: queue that stores methylation data for further processing 
    # methylation_read_number: counter for tracking number of reads processed
    # bams_analysed: counter for tracking number of BAM files processed
    # runs: dictionary for run IDs and dates from BAM header 

    while True:
        # retrieve next BAM file from queue 
        bamfile = file_queue.get()
        if bamfile is None:
            break
        # initialize defaultdict for alignment data 
        alignm = defaultdict(list) 

        # generate a BAM index if not already there
        if not os.path.exists(bamfile + '.bai'):
            pysam.index(bamfile)

        # open BAM file 
        bam = pysam.AlignmentFile(bamfile)

        # get run ID and date from first read group in BAM file 
        runs[bam.header.as_dict()['RG'][0]['ID']] = bam.header.as_dict()['RG'][0]['DT']

        # iterating over each alignment in BAM file
        reads_seen = 0
        for alignment in bam:
            reads_seen += 1
            # ignoring secondary and supplementary alignments 
            if not alignment.is_secondary and not alignment.is_supplementary:
                if alignment.mapping_quality >= 10:
                    # collect relevant alignment data 
                    alignm[alignment.qname].append(
                        (
                            alignment.reference_name,
                            alignment.get_aligned_pairs(),
                            alignment.modified_bases,
                            alignment.is_reverse,
                            alignment.get_tag('st'),
                            alignment.get_tag('RG'),
                            alignment.get_tag('rn'),
                            alignment.qname,
                            alignment.get_tag('qs'),
                            alignment.infer_read_length(),
                            alignment.mapping_quality
                        )
                    )
            if reads_seen % read_amount == 0 and alignm:
                unique_alignments = list(iter_unique_methylation_alignments(alignm))
                for chunk in chunk_records(unique_alignments, chunk_size):
                    methylation_queue.put(chunk)
                    with methylation_read_number.get_lock():
                        methylation_read_number.value += len(chunk)
                alignm.clear()

        if alignm:
            unique_alignments = list(iter_unique_methylation_alignments(alignm))
            for chunk in chunk_records(unique_alignments, chunk_size):
                methylation_queue.put(chunk)
                with methylation_read_number.get_lock():
                    methylation_read_number.value += len(chunk)
        bam.close()
        with bams_analysed.get_lock():
            bams_analysed.value += 1

def process_methylation_data(read_row, sites_index):
    # create alignment pairs panda df
    aligned_pairs = pd.DataFrame(read_row[1], columns=["query_pos", "ref_pos"])
    modified_bases = read_row[2]

    chrom_sites = sites_index.get(read_row[0])
    if chrom_sites is None:
        return []

    if read_row[3] is True:  # if read is reversed
        modified_positions = modified_bases.get(("C", 1, "m"))
        if not modified_positions:
            return []
        methylation_pos_df = pd.DataFrame(
            modified_positions,
            columns=["query_pos", "methylation"],
        )
        res = (
            aligned_pairs.join(methylation_pos_df.set_index("query_pos"), on="query_pos")
            .dropna()
            .reset_index(drop=True)
        )
        cpgs = chrom_sites["reverse"].join(res.set_index("ref_pos"), how="inner")[
            ["epic_id", "methylation"]
        ]
    else:
        modified_positions = modified_bases.get(("C", 0, "m"))
        if not modified_positions:
            return []
        methylation_pos_df = pd.DataFrame(
            modified_positions,
            columns=["query_pos", "methylation"],
        )
        res = (
            aligned_pairs.join(methylation_pos_df.set_index("query_pos"), on="query_pos")
            .dropna()
            .reset_index(drop=True)
        )
        cpgs = chrom_sites["forward"].join(res.set_index("ref_pos"), how="inner")[
            ["epic_id", "methylation"]
        ]

    if cpgs.empty:
        return []

    cpgs["methylation"] = cpgs["methylation"] / 255
    cpgs["scores_per_read"] = len(cpgs.index)
    cpgs["binary_methylation"] = (cpgs["methylation"] >= 0.8).astype(int)
    cpgs["read_id"] = read_row[7]
    cpgs["start_time"] = read_row[4]
    cpgs["run_id"] = read_row[5]
    cpgs["QS"] = read_row[8]
    cpgs["read_length"] = read_row[9]
    cpgs["map_qs"] = read_row[10]
    return cpgs[
        [
            "epic_id",
            "methylation",
            "scores_per_read",
            "binary_methylation",
            "read_id",
            "start_time",
            "run_id",
            "QS",
            "read_length",
            "map_qs",
        ]
    ].values.tolist()


def methylation_reader(
    methylation_queue,
    result_queue,
    sites_index,
    finished,
    methylation_read_number_analysed,
):
    # function for parsing methylation information from the sequencing data 

    # params: 
    # methylation_queue: queue with chunks of methylation data for further processing 
    # methylation_list: a list for storing processed methylation data 
    # sites_set: pd dataframe with CpG site annotation for sites of interest 
    # finished: counter for completion of parallel tasks 
    # methylation_read_number_analysed: counter tracking number of processed reads 

    while True:
        # get the next chunk of methylation data from queue 
        data = methylation_queue.get()

        if data is None:
            time.sleep(1)
            with finished.get_lock():
                finished.value += 1
            break

        rows = []
        for read_row in data:
            with methylation_read_number_analysed.get_lock():
                methylation_read_number_analysed.value += 1
            rows.extend(process_methylation_data(read_row, sites_index))
        if rows:
            result_queue.put(rows)
    result_queue.put(None)

def progress(finished,methylation_threads,bams_analysed,bams_amount,methylation_read_number_analysed,methylation_read_number):
    # function to track progress of parallel tasks 
    while True:
        if finished.value >= methylation_threads:
            break
        time.sleep(1)
        print('\r' + 'Bams: ' + str(bams_analysed.value) + '/' + str(bams_amount.value) + ' Reads: ' + str(methylation_read_number_analysed.value) + '/' + str(methylation_read_number.value), end="")


def main(inputs, recursive, io_threads, methylation_threads, sites, sample, output, bam_filter, read_amount):
    # 

    # params:
    # inputs: directory where BAM files are looked for  
    # recusive: boolean whether searching for files should be recursive
    # io_threads: number of threads used for io-handling
    # methylation_threads: number of threads used for methylation data extraction 
    # sites: BED file with CpG site annotation for sites of interest 
    # sample: sample file name 
    # output: output directory 
    # bam_filter: optional string for filtering on barcodes 

    # initializing variables for multiprocessing 
    file_queue = Queue()
    manager = Manager()
    methylation_queue = Queue(100)
    result_queue = Queue(100)
    finished = Value('i',0,lock=True)
    methylation_read_number = Value('i',0,lock=True)
    methylation_read_number_analysed = Value('i',0,lock=True)
    bams_amount = Value('i',0,lock=True)
    bams_analysed = Value('i',0,lock=True)
    runs = manager.dict()

    # loading CpG sites of interest with annotion from BED file, convert to panda df 
    print('Loading Sites set.')
    sites_set = pd.read_csv(sites , sep='\t', index_col=False,
            names=["chromosome", "start", "end", "epic_id"],
            dtype={"chromosome": str, "start": np.int32, "end": np.int32, "epic_id": str},
    )
    sites_set['end'] = sites_set['end']-1
    sites_index = build_sites_index(sites_set)
    # starting to process the reading of BAM files and enqueue their paths
    print('Start getting files.')
    bam_process = Process(target=read_bam_files, args=(inputs, recursive, file_queue, bams_amount, bam_filter))
    bam_process.start()
    print('Start reading files.')
    print()
    # initialize and start IO handler for processing of BAM files 
    io_processes = []
    chunk_size = max(1, read_amount // max(1, methylation_threads))
    for i in range(io_threads):
        io_processes.append(Process(target=io_handler, args=(file_queue, methylation_queue, methylation_read_number, bams_analysed, runs, read_amount, chunk_size)))
        io_processes[i].start()
    # initialize and start methylation reader for processing of methylation data 
    methylation_processes = []
    for i in range(methylation_threads):
        methylation_processes.append(Process(target=methylation_reader, args=(methylation_queue, result_queue, sites_index, finished, methylation_read_number_analysed)))
        methylation_processes[i].start()

    # track progress of analysis 
    progress_process = Process(target=progress,args=(finished,methylation_threads,bams_analysed,bams_amount,methylation_read_number_analysed,methylation_read_number))
    progress_process.start()

    # waiting for BAM processing to complete
    bam_process.join()

    for i in range(io_threads):
        file_queue.put(None)

    # waiting for IO handling processes to complete
    for i in range(io_threads):
        io_processes[i].join()

    for j in range(methylation_threads):
        methylation_queue.put(None)

    # waiting for methylataion reader processes to complete
    methylation_rows = []
    finished_readers = 0
    while finished_readers < methylation_threads:
        result = result_queue.get()
        if result is None:
            finished_readers += 1
        else:
            methylation_rows.extend(result)

    for i in range(methylation_threads):
        methylation_processes[i].join()

    progress_process.join()

    # saving methylation data, if any was collected 
    print()
    if methylation_rows:
        print('Saving data to feather.')
        # creating panda df from methylation list 
        methylation = pd.DataFrame(methylation_rows, columns=["epic_id", "methylation", "scores_per_read", "binary_methylation", "read_id", "start_time", "run_id","QS","read_length","map_qs"]).sort_values('start_time')

        # for each run, adjust start times based on run-specific metadata 
        methylation_runs = []
        for run_id, st in runs.items():

            t1 = arrow.get(methylation[methylation['run_id']==run_id].sort_values('start_time').iloc[0]['start_time'])
            tdif = np.floor((t1-arrow.get(st)).total_seconds()/3600)*3600
            run = methylation[methylation['run_id']==run_id]
            run.loc[:,'start_time'] = run['start_time'].apply(lambda a: int((arrow.get(a)-arrow.get(st)).total_seconds()-tdif))
            methylation_runs.append(run)

        methylation = pd.concat(methylation_runs)
        # ensure output dir exists 
        if not os.path.exists(output):
            os.makedirs(output)
        # save the methylation data to a Feather file for efficient storage and access
        methylation.sort_values('start_time').reset_index(drop=True).to_feather(output + '/' + sample + '.feather')
    else:
        print('No data retrieved.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # command line arguments that scripts accepts 
    parser.add_argument('-i', '--inputs', type=str, default='.', required=True, help='Filepath of BAM files')
    parser.add_argument('-r', '--recursive', default=False, action='store_true', help='Recursively monitor subdirectories')
    parser.add_argument('--io_threads', type=int, default=2, help='Number of Threads used for io-handling')
    parser.add_argument('--methylation_threads', default=4, type=int, help='Number of Threads used for methylation data extraction')
    parser.add_argument('--sites', type=str, default='.', required=True, help='File with CpG-Sites-Annotation in bed format')
    parser.add_argument('-s', '--sample', type=str, required=True, help='Name of the Sample')
    parser.add_argument('-o', '--output', type=str,required=True, help="Path to output foldere.")
    parser.add_argument('--filter', type=str, default=None, help='String that needs to be present in path. Useful for filtering on barcodes.')
    parser.add_argument('--read_amount', type=int, default=4000, help='Number of reads to process per batch.')
    args = parser.parse_args()
    # executing main function with parsed arguments 
    main(args.inputs, args.recursive, args.io_threads, args.methylation_threads, args.sites, args.sample, args.output, args.filter, args.read_amount)
