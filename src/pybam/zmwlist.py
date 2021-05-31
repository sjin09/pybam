import gzip
import pysam
import logging
from pybam.bamlib import BAM


def bam2zmw(infile, outfile):
    zmw_lst = []
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        zmw_lst.append(read.zmw)
    outfile.write("{}\n".format("\n".join(zmw_lst)))


def fq2zmw(infile, outfile):
    zmw_lst = []
    seqfile = open(infile) if infile.endswith((".fq", ".fastq")) else gzip.open(infile)
    for i, j in enumerate(seqfile):
        k = i % 4
        if k == 0:  # header
            seq_id = j.strip()
            seq_id = seq_id if not isinstance(seq_id, bytes) else seq_id.decode("utf-8")
            seq_id = seq_id.replace("@", "")
            zmw = seq_id.split("/")[1]
            zmw_lst.append(zmw)
        elif k == 1:  # sequence
            continue
        elif k == 2:
            continue  # plus
        elif k == 3:  # quality
            continue
    outfile.write("{}\n".format("\n".join(zmw_lst)))


def zmwlist(infile, outfile):
    if infile.endswith(".bam"):
        bam2zmw(infile, outfile)
    elif infile.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
        fq2zmw(infile, outfile)
    else:
        logging.error("pybam does not support the provided INPUT file")
