import sys
import json
import gzip
import pysam
import logging
import natsort
import numpy as np
from pybam.bamClass import BAM


def bam_stat(bamfile):
    ccs_hsh = {}
    alignment_file = pysam.AlignmentFile(bamfile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        bq_sum = sum(read.bq_int_lst)
        qv = bq_sum / float(read.qlen)
        hq_count = read.bq_ascii.count("~")
        hq_proportion = "{:.2f}".format((hq_count * 100) / float(read.qlen))
        ccs_hsh[read.zmw] = [qv, hq_proportion, read.qlen]


def fq_stat(fqfile):
    ccs_hsh = {}
    seqfile = open(fqfile) if fqfile.endswith((".fq", ".fastq")) else gzip.open(fqfile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  # header
            seq_id = j.strip()
            seq_id = str(seq_id) if not isinstance(seq_id, bytes) else seq_id.decode("utf-8") 
            seq_id = seq_id.replace("@", "")
            zmw = seq_id.split("/")[1]
        elif k == 1: # sequence
            seq = j.strip()
            seq = str(seq) if not isinstance(seq, bytes) else seq.decode("utf-8") 
            seq_len = len(seq)
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            seq_bq = j.strip()
            seq_bq = str(seq_bq) if not isinstance(seq_bq, bytes) else seq_bq.decode("utf-8") 
            seq_bq_lst = [ord(_bq) - 33 for _bq in list(seq_bq)]
            bq_sum = sum(seq_bq_lst)
            qv = bq_sum/float(seq_len)
            hq_count = seq_bq.count("~")
            hq_proportion = "{:.2f}".format((hq_count * 100)/float(seq_len)) 
            ccs_hsh[zmw] = [qv, hq_proportion, seq_len]
    return ccs_hsh


def ccs_stat(ccs):
    if ccs.endswith(".bam"):
        ccs_hsh = bam_stat(ccs)
    elif ccs.ensdwith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        ccs_hsh = fq_stat(ccs)
    return ccs_hsh


def subread_stat(subread):
    state = 0
    subread_hsh = {}
    alignment_file = pysam.AlignmentFile(subread, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        if state == 0:
            state = 1
            zmw = read.zmw
            length_lst = [read.qlen]
        else:
            nzmw = read.zmw
            if zmw == nzmw:
                length_lst.append(read.qlen)
            else:
                subread_hsh[zmw] = length_lst
                state = 0
                zmw = read.zmw
                length_lst = [read.qlen]
    return subread_hsh


def return_zmwstat(ccs_hsh, subread_hsh, outfile):
    # return: header
    hstr = "ZMW QV Q93(%) CCS_LENGTH NORMAL_COUNT FRAGMENT_COUNT CHIMIERA_COUNT SUBREAD_COUNT \
        MIN_SUBREAD_LENGTH MEDIAN_SUBREAD_LENGTH MAX_SUBREAD_LENGTH \
        LOWER_SUBREAD_THRESHOLD UPPER_SUBREAD_THRESHOLD SUBREAD_LENGTHS"
    hstr = "\t".join(hstr.split())
    outfile.write("{}\n".format(hstr))

    # return: zmw statistics
    zmw_lst = natsort.natsorted(list(subread_hsh.keys()))
    for zmw in zmw_lst:
        # ccs
        ccs_qv = "."
        ccs_length = "."
        ccs_hq_proportion = "."
        if zmw in ccs_hsh:
            ccs_qv, ccs_hq_proportion, ccs_length = ccs_hsh[zmw]
        # subreads
        subread_length_lst = subread_hsh[zmw]
        subread_count = len(subread_length_lst)
        min_subread_length = min(subread_length_lst)
        median_subread_length = np.median(subread_length_lst)
        max_subread_length = max(subread_length_lst)
        lower_subread_threshold = 0.5 * median_subread_length
        upper_subread_threshold = 2 * median_subread_length
        normal_subread_hsh = {
            i: j
            for i, j in enumerate(subread_length_lst)
            if j > lower_subread_threshold and j < upper_subread_threshold
        }
        fragmented_subread_hsh = {
            i: j
            for i, j in enumerate(subread_length_lst)
            if j < lower_subread_threshold
        }
        chimera_subread_hsh = {
            i: j
            for i, j in enumerate(subread_length_lst)
            if j > upper_subread_threshold
        }
        normal_count = len(normal_subread_hsh)
        chimera_count = len(chimera_subread_hsh)
        fragmented_count = len(fragmented_subread_hsh)
        subread_length_lst = [str(_length) for _length in subread_length_lst]
        subread_lengths = ",".join(subread_length_lst)

        ## return:
        outstr = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            zmw,
            ccs_qv,
            ccs_hq_proportion,
            ccs_length,
            normal_count,
            fragmented_count,
            chimera_count,
            subread_count,
            min_subread_length,
            median_subread_length,
            max_subread_length,
            lower_subread_threshold,
            upper_subread_threshold,
            subread_lengths,
        )
        outfile.write("{}".format(outstr))


def zmwstat(ccs, subreads, outfile):
    if ccs.endswith(".bam") and subreads.endswith(".bam"):
        ccs_hsh = ccs_stat(ccs)
        subread_hsh = subread_stat(subreads)
        return_zmwstat(ccs_hsh, subread_hsh, outfile)
    elif ccs.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")) and subreads.endswith(
        ".bam"
    ):
        ccs_hsh = ccs_stat(ccs)
        subread_hsh = subread_stat(subreads)
        return_zmwstat(ccs_hsh, subread_hsh, outfile)

    else:
        logging.error("zmwstat doesn't support the provided input files")
