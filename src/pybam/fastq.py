import os
import gzip
import pysam
import logging
from pybam.bamClass import BAM


def load_blacklist(blacklist):
    blacklist_lst = [line.strip() for line in open(blacklist).readlines()]
    return blacklist_lst


def load_whitelist(infile, blacklist):
    blacklist_set = set(blacklist)
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    zmw_lst = set([line.query_name for line in alignment_file])
    whitelist_set = zmw_lst.difference(blacklist_set) 
    return whitelist_set


def read_indexed_bamfile(infile):
    read_indexed = pysam.IndexedReads(infile)
    read_indexed.build()
    return read_indexed


def return_fastq(infile, blacklist, outfile):
    if os.path.exists(blacklist):
        zmw_blacklist = load_blacklist(blacklist)
        state = 0 if len(zmw_blacklist) == 0 else 1
    else:
        state = 0
    print(blacklist, state)
    if state:
        fqfile = open(outfile, "w")
        alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
        for line in alignment_file:
            read = BAM(line)
            fqfile.write("@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii))
    else:
        counter = 0 
        zmw_blacklist = set(zmw_blacklist)
        zmw_indexed = read_indexed_bamfile(infile)
        zmw_whitelist = load_whitelist(infile, zmw_blacklist)
        alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
        for zmw in zmw_whitelist:
            try:
                zmw_indexed.find(zmw)
            except KeyError:
                pass
            else:
                iterator = zmw_indexed.find(zmw)
                for x in iterator:
                    print(x)
            counter += 1
            if counter == 10:
                break

def return_gzip_fastq(infile, blacklist, outfile):
    fqfile = gzip.open(outfile, "wb")
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        fqstr = "@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii)
        fqstr = fqstr.encode("utf-8")
        fqfile.write(fqstr)


# def return_whitelist_fastq(infile, blacklist, outfile):


## def return_whitelist_gzip_fastq(infile, blacklist, outfile):

def bam2fastq(infile, blacklist, outfile):
    if infile.endswith((".sam", ".bam")):
        if outfile.endswith((".fq", ".fastq")):
            return_fastq(infile, blacklist, outfile)
        elif outfile.endswith((".fq.gz", ".fastq.gz")):
            return_gzip_fastq(infile, blacklist, outfile)
        else:
            logging.error("bam2fastq does not support the provided OUTPUT file format")
    else:
        logging.error("bam2fastq does not support the provided INPUT file")
