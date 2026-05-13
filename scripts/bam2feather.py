#!/usr/bin/env python3
import argparse
import os
import pathlib
import shutil
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
            "forward": dict(zip(chrom_sites["start"], chrom_sites["epic_id"])),
            "reverse": dict(zip(chrom_sites["end"], chrom_sites["epic_id"])),
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


def chunk_unique_methylation_alignments(alignments_by_read, chunk_size):
    chunk = []
    for read_alignment in iter_unique_methylation_alignments(alignments_by_read):
        chunk.append(read_alignment)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def discover_bam_files(path, recursive, bam_filter):
    d = pathlib.Path(path)
    globber = d.rglob if recursive else d.glob
    files = []
    for file in natsorted(globber("*.bam")):
        if bam_filter is None or bam_filter in str(file):
            files.append(str(file))
    return files


def read_bam_files(path, recursive, file_queue, bams_amount, bam_filter):
    # function to read BAM files from minimap2 output  

    # params: 
    # path: directory where BAM files are looked for 
    # recursive: boolean whether searching for files should be recursive
    # file_queue: queue that stores BAM files for further processing
    # bams_amount: counter for tracking total number of found BAM files 
    # bam_filter: optional string for filtering on barcodes 

    files = discover_bam_files(path, recursive, bam_filter)
    with bams_amount.get_lock():
        bams_amount.value = len(files)
    for file in files:
        file_queue.put(file)


def ensure_bam_index(bamfile):
    if not os.path.exists(bamfile + ".bai"):
        pysam.index(bamfile)


def _bam_region_tasks(bamfile, sites_index):
    ensure_bam_index(bamfile)
    with pysam.AlignmentFile(bamfile) as bam:
        references = set(bam.references)
    regions = [chromosome for chromosome in sites_index if chromosome in references]
    if not regions:
        return [(bamfile, None)]
    return [(bamfile, chromosome) for chromosome in regions]


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
                for chunk in chunk_unique_methylation_alignments(alignm, chunk_size):
                    methylation_queue.put(chunk)
                    with methylation_read_number.get_lock():
                        methylation_read_number.value += len(chunk)
                alignm.clear()

        if alignm:
            for chunk in chunk_unique_methylation_alignments(alignm, chunk_size):
                methylation_queue.put(chunk)
                with methylation_read_number.get_lock():
                    methylation_read_number.value += len(chunk)
        bam.close()
        with bams_analysed.get_lock():
            bams_analysed.value += 1

def process_methylation_data(read_row, sites_index):
    modified_bases = read_row[2]

    chrom_sites = sites_index.get(read_row[0])
    if chrom_sites is None:
        return []

    query_to_ref = {}
    for query_pos, ref_pos in read_row[1]:
        if query_pos is None or ref_pos is None:
            continue
        query_to_ref[query_pos] = ref_pos

    if read_row[3] is True:  # if read is reversed
        modified_positions = modified_bases.get(("C", 1, "m"))
        if not modified_positions:
            return []
        epic_lookup = chrom_sites["reverse"]
    else:
        modified_positions = modified_bases.get(("C", 0, "m"))
        if not modified_positions:
            return []
        epic_lookup = chrom_sites["forward"]

    cpgs = []
    for query_pos, methylation_raw in modified_positions:
        ref_pos = query_to_ref.get(query_pos)
        if ref_pos is None:
            continue
        epic_id = epic_lookup.get(ref_pos)
        if epic_id is None:
            continue
        methylation = methylation_raw / 255
        cpgs.append([epic_id, methylation, int(methylation >= 0.8)])

    if not cpgs:
        return []

    scores_per_read = len(cpgs)
    rows = []
    for epic_id, methylation, binary_methylation in cpgs:
        rows.append(
            [
                epic_id,
                methylation,
                scores_per_read,
                binary_methylation,
                read_row[7],
                read_row[4],
                read_row[5],
                read_row[8],
                read_row[9],
                read_row[10],
            ]
        )
    return rows


def modified_c_positions(modified_bases, is_reverse):
    if is_reverse:
        return modified_bases.get(("C", 1, "m")) or []
    return modified_bases.get(("C", 0, "m")) or []


def reference_positions_for_query_positions(alignment, query_positions):
    if not query_positions or alignment.reference_start is None:
        return {}

    targets = sorted(set(query_positions))
    target_index = 0
    query_pos = 0
    reference_pos = alignment.reference_start
    mapped = {}

    for operation, length in alignment.cigartuples or []:
        if target_index >= len(targets):
            break

        if operation in (0, 7, 8):
            query_end = query_pos + length
            while target_index < len(targets) and targets[target_index] < query_end:
                target = targets[target_index]
                if target >= query_pos:
                    mapped[target] = reference_pos + (target - query_pos)
                target_index += 1
            query_pos = query_end
            reference_pos += length
        elif operation in (1, 4):
            query_end = query_pos + length
            while target_index < len(targets) and targets[target_index] < query_end:
                target_index += 1
            query_pos = query_end
        elif operation in (2, 3):
            reference_pos += length
        elif operation in (5, 6):
            continue

    return mapped


def process_alignment_methylation_data(alignment, sites_index):
    chrom_sites = sites_index.get(alignment.reference_name)
    if chrom_sites is None:
        return []

    modified_positions = modified_c_positions(alignment.modified_bases, alignment.is_reverse)
    if not modified_positions:
        return []

    epic_lookup = chrom_sites["reverse"] if alignment.is_reverse else chrom_sites["forward"]
    query_positions = [query_pos for query_pos, _ in modified_positions]
    query_to_ref = reference_positions_for_query_positions(alignment, query_positions)

    cpgs = []
    for query_pos, methylation_raw in modified_positions:
        ref_pos = query_to_ref.get(query_pos)
        if ref_pos is None:
            continue
        epic_id = epic_lookup.get(ref_pos)
        if epic_id is None:
            continue
        methylation = methylation_raw / 255
        cpgs.append([epic_id, methylation, int(methylation >= 0.8)])

    if not cpgs:
        return []

    scores_per_read = len(cpgs)
    start_time = alignment.get_tag("st")
    run_id = alignment.get_tag("RG")
    qs = alignment.get_tag("qs")
    read_length = alignment.infer_read_length()
    mapping_quality = alignment.mapping_quality
    return [
        [
            epic_id,
            methylation,
            scores_per_read,
            binary_methylation,
            alignment.qname,
            start_time,
            run_id,
            qs,
            read_length,
            mapping_quality,
        ]
        for epic_id, methylation, binary_methylation in cpgs
    ]


def iter_primary_methylation_alignments(bam, region=None):
    source = bam.fetch(region) if region is not None else bam
    for alignment in source:
        if alignment.is_secondary or alignment.is_supplementary:
            continue
        if alignment.mapping_quality < 10:
            continue
        yield alignment


def write_pending_rows(chunk_dir, worker_label, chunk_index, pending_rows):
    if not pending_rows:
        return chunk_index
    write_chunk_file(chunk_dir, worker_label, chunk_index, pending_rows)
    return chunk_index + 1


def process_alignment_batch(
    alignments_by_read,
    sites_index,
    pending_rows,
    chunk_dir,
    chunk_row_target,
    worker_label,
    chunk_index,
    methylation_read_number,
    methylation_read_number_analysed,
):
    for alignments in alignments_by_read.values():
        if len(alignments) != 1:
            continue
        alignment = alignments[0]
        if alignment.modified_bases == {}:
            continue
        with methylation_read_number.get_lock():
            methylation_read_number.value += 1
        rows = process_alignment_methylation_data(alignment, sites_index)
        with methylation_read_number_analysed.get_lock():
            methylation_read_number_analysed.value += 1
        pending_rows.extend(rows)
        if len(pending_rows) >= chunk_row_target:
            chunk_index = write_pending_rows(chunk_dir, worker_label, chunk_index, pending_rows)
            pending_rows = []
    alignments_by_read.clear()
    return pending_rows, chunk_index


def methylation_task_worker(
    task_queue,
    sites_index,
    finished,
    methylation_read_number,
    methylation_read_number_analysed,
    tasks_analysed,
    runs,
    chunk_dir,
    chunk_row_target,
    read_amount,
    worker_label,
):
    pending_rows = []
    chunk_index = 0

    while True:
        task = task_queue.get()
        if task is None:
            chunk_index = write_pending_rows(chunk_dir, worker_label, chunk_index, pending_rows)
            pending_rows = []
            with finished.get_lock():
                finished.value += 1
            break

        bamfile, region = task
        with pysam.AlignmentFile(bamfile) as bam:
            for read_group in bam.header.as_dict().get("RG", []):
                if "ID" in read_group and "DT" in read_group:
                    runs[read_group["ID"]] = read_group["DT"]

            alignments_by_read = defaultdict(list)
            reads_seen = 0
            read_batch_size = max(1, read_amount)
            for alignment in iter_primary_methylation_alignments(bam, region):
                alignments_by_read[alignment.qname].append(alignment)
                reads_seen += 1
                if reads_seen % read_batch_size == 0 and alignments_by_read:
                    pending_rows, chunk_index = process_alignment_batch(
                        alignments_by_read,
                        sites_index,
                        pending_rows,
                        chunk_dir,
                        chunk_row_target,
                        worker_label,
                        chunk_index,
                        methylation_read_number,
                        methylation_read_number_analysed,
                    )

            if alignments_by_read:
                pending_rows, chunk_index = process_alignment_batch(
                    alignments_by_read,
                    sites_index,
                    pending_rows,
                    chunk_dir,
                    chunk_row_target,
                    worker_label,
                    chunk_index,
                    methylation_read_number,
                    methylation_read_number_analysed,
                )

        with tasks_analysed.get_lock():
            tasks_analysed.value += 1


def methylation_reader(
    methylation_queue,
    sites_index,
    finished,
    methylation_read_number_analysed,
    chunk_dir,
    chunk_row_target,
    worker_label,
):
    # function for parsing methylation information from the sequencing data 

    # params: 
    # methylation_queue: queue with chunks of methylation data for further processing 
    # methylation_list: a list for storing processed methylation data 
    # sites_set: pd dataframe with CpG site annotation for sites of interest 
    # finished: counter for completion of parallel tasks 
    # methylation_read_number_analysed: counter tracking number of processed reads 

    pending_rows = []
    chunk_index = 0

    while True:
        # get the next chunk of methylation data from queue 
        data = methylation_queue.get()

        if data is None:
            if pending_rows:
                write_chunk_file(chunk_dir, worker_label, chunk_index, pending_rows)
            with finished.get_lock():
                finished.value += 1
            break

        for read_row in data:
            with methylation_read_number_analysed.get_lock():
                methylation_read_number_analysed.value += 1
            pending_rows.extend(process_methylation_data(read_row, sites_index))
            if len(pending_rows) >= chunk_row_target:
                write_chunk_file(chunk_dir, worker_label, chunk_index, pending_rows)
                chunk_index += 1
                pending_rows = []


RESULT_COLUMNS = [
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


def write_chunk_file(chunk_dir, worker_label, chunk_index, rows):
    os.makedirs(chunk_dir, exist_ok=True)
    chunk_path = os.path.join(chunk_dir, f"{worker_label}-chunk-{chunk_index:05d}.feather")
    pd.DataFrame(rows, columns=RESULT_COLUMNS).to_feather(chunk_path)
    return chunk_path


def normalize_run_start_times(methylation, runs):
    methylation_runs = []
    for run_id, st in runs.items():
        run = methylation[methylation['run_id'] == run_id].copy()
        if run.empty:
            continue

        t1 = arrow.get(run.sort_values('start_time').iloc[0]['start_time'])
        tdif = np.floor((t1 - arrow.get(st)).total_seconds() / 3600) * 3600
        normalized = run['start_time'].apply(
            lambda a: int((arrow.get(a) - arrow.get(st)).total_seconds() - tdif)
        )
        run = run.assign(start_time=normalized.astype('int64'))
        methylation_runs.append(run)

    if not methylation_runs:
        return methylation

    return pd.concat(methylation_runs, ignore_index=True)


def progress(finished, worker_count, tasks_analysed, tasks_amount, methylation_read_number_analysed, methylation_read_number):
    # function to track progress of parallel tasks 
    started = time.time()
    last_newline = started
    while True:
        if finished.value >= worker_count:
            break
        time.sleep(1)
        elapsed = max(time.time() - started, 1e-6)
        analysed = methylation_read_number_analysed.value
        queued = methylation_read_number.value
        message = (
            f"Tasks: {tasks_analysed.value}/{tasks_amount.value} "
            f"Reads: {analysed}/{queued} "
            f"Rate: {analysed / elapsed:.1f} reads/s "
            f"Elapsed: {elapsed:.0f}s"
        )
        if time.time() - last_newline >= 30:
            print(message, flush=True)
            last_newline = time.time()
        else:
            print('\r' + message, end="", flush=True)


def main(inputs, recursive, io_threads, methylation_threads, sites, sample, output, bam_filter, read_amount):
    # 

    # params:
    # inputs: directory where BAM files are looked for  
    # recusive: boolean whether searching for files should be recursive
    # io_threads: retained for command-line compatibility
    # methylation_threads: number of BAM/contig workers used for methylation data extraction
    # sites: BED file with CpG site annotation for sites of interest 
    # sample: sample file name 
    # output: output directory 
    # bam_filter: optional string for filtering on barcodes 

    # initializing variables for multiprocessing 
    task_queue = Queue()
    manager = Manager()
    finished = Value('i',0,lock=True)
    methylation_read_number = Value('i',0,lock=True)
    methylation_read_number_analysed = Value('i',0,lock=True)
    tasks_amount = Value('i',0,lock=True)
    tasks_analysed = Value('i',0,lock=True)
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
    bam_files = discover_bam_files(inputs, recursive, bam_filter)
    tasks = []
    for bam_file in bam_files:
        tasks.extend(_bam_region_tasks(bam_file, sites_index))
    with tasks_amount.get_lock():
        tasks_amount.value = len(tasks)
    for task in tasks:
        task_queue.put(task)
    print('Start reading files.')
    print()
    chunk_size = max(1, read_amount // max(1, methylation_threads))
    chunk_dir = os.path.join(output, ".bam2feather-chunks")
    chunk_row_target = max(100000, chunk_size * 200)
    worker_count = max(1, methylation_threads)
    methylation_processes = []
    for i in range(worker_count):
        methylation_processes.append(
            Process(
                target=methylation_task_worker,
                args=(
                    task_queue,
                    sites_index,
                    finished,
                    methylation_read_number,
                    methylation_read_number_analysed,
                    tasks_analysed,
                    runs,
                    chunk_dir,
                    chunk_row_target,
                    read_amount,
                    f"worker-{i:02d}",
                ),
            )
        )
        methylation_processes[i].start()

    # track progress of analysis 
    progress_process = Process(target=progress,args=(finished,worker_count,tasks_analysed,tasks_amount,methylation_read_number_analysed,methylation_read_number))
    progress_process.start()

    for _ in range(worker_count):
        task_queue.put(None)

    for i in range(worker_count):
        methylation_processes[i].join()

    progress_process.join()

    # saving methylation data, if any was collected 
    print()
    chunk_paths = sorted(pathlib.Path(chunk_dir).glob("*.feather"))
    if chunk_paths:
        print('Saving data to feather.')
        methylation = pd.concat(
            (pd.read_feather(chunk_path) for chunk_path in chunk_paths),
            ignore_index=True,
        ).sort_values('start_time')

        # for each run, adjust start times based on run-specific metadata 
        methylation = normalize_run_start_times(methylation, runs)
        # ensure output dir exists 
        if not os.path.exists(output):
            os.makedirs(output)
        # save the methylation data to a Feather file for efficient storage and access
        methylation.sort_values('start_time').reset_index(drop=True).to_feather(output + '/' + sample + '.feather')
        shutil.rmtree(chunk_dir, ignore_errors=True)
    else:
        print('No data retrieved.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # command line arguments that scripts accepts 
    parser.add_argument('-i', '--inputs', type=str, default='.', required=True, help='Filepath of BAM files')
    parser.add_argument('-r', '--recursive', default=False, action='store_true', help='Recursively monitor subdirectories')
    parser.add_argument('--io_threads', type=int, default=2, help='Retained for command-line compatibility with earlier releases')
    parser.add_argument('--methylation_threads', default=4, type=int, help='Number of BAM/contig workers used for methylation data extraction')
    parser.add_argument('--sites', type=str, default='.', required=True, help='File with CpG-Sites-Annotation in bed format')
    parser.add_argument('-s', '--sample', type=str, required=True, help='Name of the Sample')
    parser.add_argument('-o', '--output', type=str,required=True, help="Path to output foldere.")
    parser.add_argument('--filter', type=str, default=None, help='String that needs to be present in path. Useful for filtering on barcodes.')
    parser.add_argument('--read_amount', type=int, default=4000, help='Number of reads to process per batch.')
    args = parser.parse_args()
    # executing main function with parsed arguments 
    main(args.inputs, args.recursive, args.io_threads, args.methylation_threads, args.sites, args.sample, args.output, args.filter, args.read_amount)
