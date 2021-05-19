import os
import gzip
import pysam
import logging
from pybam.bamClass import BAM


def return_fasta(infile, outfile):
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        outfile.write("@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii))

# def return_gzip_fasta(infile, outfile):
#     fafile = gzip.open(outfile, "wb")
#     alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
#     for line in alignment_file:
#         read = BAM(line)
#         fastr = "@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii)
#         fastr = fastr.encode("utf-8")
#         fafile.write(fastr)

def bam2fasta(infile, outfile):
    if infile.endswith((".sam", ".bam")):
        return_fasta(infile, outfile)
    else:
        logging.error("bam2fastq does not support the provided INPUT file")
