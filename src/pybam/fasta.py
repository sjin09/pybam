import os
import gzip
import pysam
import logging
from fastx.bamClass import BAM
from fastx.util import chunkstring


def return_fasta(infile, outfile):
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        outfile.write(">{}\n".format(read.qname))
        for chunk in chunkstring(read.seq):
            outfile.write("{}\n".format(chunk))


def bam2fasta(infile, outfile):
    if infile.endswith((".sam", ".bam")):
        return_fasta(infile, outfile)
    else:
        logging.error("bam2fastq does not support the provided INPUT file")

