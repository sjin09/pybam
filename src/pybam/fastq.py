import gzip
import pysam
import logging
from pybam.bamClass import BAM


def return_fastq(infile, outfile):
    fqfile = open(outfile, "w")
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        fqfile.write("@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii))


def return_gzip_fastq(infile, outfile):
    fqfile = gzip.open(outfile, "wb")
    alignment_file = pysam.AlignmentFile(infile, "rb", check_sq=False)
    for line in alignment_file:
        read = BAM(line)
        fqstr = "@{}\n{}\n+\n{}\n".format(read.qname, read.qseq, read.bq_ascii)
        fqstr = fqstr.encode("utf-8")
        fqfile.write(fqstr)


def bam2fastq(infile, outfile):
    if infile.endswith((".sam", ".bam")):
        if outfile.endswith((".fq", ".fastq")):
            return_fastq(infile, outfile)
        elif outfile.endswith((".fq.gz", ".fastq.gz")):
            return_gzip_fastq(infile, outfile)
        else:
            logging.error("pybam does not support the provided INPUT file")
    else:
        logging.error("pybam does not support the provided INPUT file")
