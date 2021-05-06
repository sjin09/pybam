from os import kill
import sys
import pysam
import logging
import natsort
import numpy as np
from pybam.bamClass import BAM


def ccs_stat():
    ccs_hash = {}
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  # header
            seq_id = j.strip()
            seq_id = seq_id if not isinstance(seq_id, bytes) else seq_id.decode("utf-8")
            seq_id = seq_id.replace("@", "")
            zmw = seq_id.split("/")[1]
        elif k == 1:  # sequence
            seq = j.strip()
            seq = seq if not isinstance(seq, bytes) else seq.decode("utf-8")
            seq_len = len(seq)
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            bq_ascii = j.strip()
            bq_ascii = bq_ascii if not isinstance(bq_ascii, bytes) else bq_ascii.decode("utf-8")
            hq_base_count = bq_ascii.count("~")
            hq_base_fraction = "{:.2f}".format(hq_base_count/float(seq_len))
            ccs_hash[zmw] = [hq_base_fraction, seq_len]
    return ccs_hash


def subread_stat(subreads):
    subread_hash = {}
    status = 0
    counter = 0
    pbi = open(subreads.replace(".gz", "") + ".pbi", "w")
    sequences = (
        SeqIO.parse(subreads, "fasta")
        if subreads.endswith(".fasta")
        else SeqIO.parse(gzip.open(subreads, "rt"), "fasta")
    )
    for i in sequences:
        if status == 0:
            status = 1
            start = counter
            zmw = i.id.split("/")[1]
        else:
            next_zmw = i.id.split("/")[1]
            status = 1 if zmw == next_zmw else 0

        if status == 0:
            end = counter
            pbi.write("{}\t{}\t{}\n".format(zmw, start, end))
        counter += 1


def return_stat(ccs_hsh, subread_hsh, outfile):
    zmw_lst = list(subread_hsh.keys())
    for zmw in zmw_lst:
        ccs_bq, ccs_length = ccs_hsh[zmw]
        subread_length_lst = subread_hsh[zmw]
        subread_count = len(subread_length_lst)
        subread_lengths = ",".join(subread_length_lst)
        min_subread_length = min(subread_length_lst)
        max_subread_length = max(subread_length_lst)
        median_subread_length = np.median(subread_length_lst)
        outstr = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            zmw,
            ccs_length,
            subread_count, 
            min_subread_length,
            median_subread_length,
            max_subread_length,
            subread_lengths,
        )
        outfile.write("{}".format(outstr))


def zmwstat(ccs, subreads, outfile):
    if ccs.endswith(".bam") and subreads.endswith(".bam"):
        ccs_hsh = ccs_stat(ccs)
        # subread_hsh = subread_stat(subreads)
        # return_stat(ccs_hsh, subread_hsh, outfile)
    else:
        logging.error("zmwstat doesn't support the provided input files")
