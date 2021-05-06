from os import kill
import sys
import pysam
import logging
import natsort
import numpy as np
from pybam.bamClass import BAM


def ccs_stat(ccs):
    ccs_hsh = {}
    alignment_file = pysam.AlignmentFile(ccs, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        hq_base_count = read.bq_ascii.count("~")
        hq_base_fraction = "{:.2f}".format(hq_base_count/float(read.qlen))
        ccs_hsh[read.zmw] = [hq_base_fraction, read.qlen]
    return ccs_hsh


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
