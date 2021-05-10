import sys
import json
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
        bq_sum = sum(read.bq_int_lst)
        bq_average = bq_sum/float(read.qlen)
        hbq_count = read.bq_ascii.count("~")
        hbq_proportion = "{:.2f}".format((hbq_count * 100) / float(read.qlen))
        ccs_hsh[read.zmw] = [bq_average, hbq_proportion, read.qlen]
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


def return_stat(ccs_hsh, subread_hsh, outfile):
    # return: header
    hstr = "ZMW Q93(%) CCS_LENGTH SUBREAD_COUNT NORMAL_COUNT FRAGMENT_COUNT CHIMIERA_COUNT \
        MIN_SUBREAD_LENGTH MEDIAN_SUBREAD_LENGTH MAX_SUBREAD_LENGTH \
        LOWER_SUBREAD_THRESHOLD UPPER_SUBREAD_THRESHOLD SUBREAD_LENGTHS"
    hstr = "\t".join(hstr.split())
    outfile.write("{}\n".format(hstr))

    # return: zmw statistics
    zmw_lst = natsort.natsorted(list(subread_hsh.keys()))
    for zmw in zmw_lst:
        # ccs
        ccs_length = "."
        ccs_bq_average = "."
        ccs_hbq_proprtion = "."
        if zmw in ccs_hsh:
            ccs_bq_average, ccs_hbq_proprtion, ccs_length = ccs_hsh[zmw]

        # subreads
        subread_length_lst = subread_hsh[zmw]
        subread_count = len(subread_length_lst)
        min_subread_length = min(subread_length_lst)
        median_subread_length = np.median(subread_length_lst)
        max_subread_length = max(subread_length_lst)
        subread_length_lst = [str(_length) for _length in subread_length_lst]
        subread_lengths = ",".join(subread_length_lst)
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

        outstr = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            zmw,
            ccs_bq_average,
            ccs_hbq,
            ccs_length,
            subread_count,
            normal_count,
            fragmented_count,
            chimera_count,
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
        # print(json.dumps(ccs_hsh, indent=4))
        subread_hsh = subread_stat(subreads)
        # print(json.dumps(subread_hsh, indent=4))
        return_stat(ccs_hsh, subread_hsh, outfile)
    else:
        logging.error("zmwstat doesn't support the provided input files")
