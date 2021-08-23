#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"


# modules
import os
import sys
import logging
from pybam.fasta import bam2fasta
from pybam.fastq import bam2fastq
from pybam.zmwstat import zmwstat
from pybam.zmwlist import zmwlist
from pybam.zmw2fasta import zmw2fasta
from pybam.parse_args import parse_args

def main():
    parser, options = parse_args(program_version=__version__)

    if options.sub == "zmwlist":  # parse BAM or FASTQ file and return a list of ZMW
        zmwlist(options.input, options.output)
    elif options.sub == "zmwstat":  # generate ZMW statistisc from CCS and subreads BAM files
        zmwstat(options.ccs, options.subreads, options.output)
    elif options.sub == "zmw2fasta":  # convert BAM file to FASTQ file
        zmw2fasta(options.reads, options.dir)
    elif options.sub == "bam2fastq":  # convert BAM file to FASTQ file
        bam2fastq(options.input, options.output)
    elif options.sub == "bam2fasta":  # convert BAM file to FASTQ file
        bam2fasta(options.input, options.output)
    else:
        logging.warning("The subcommand does not exist!\n")
        parser.print_help()

if __name__ == "__main__":
    main()
    sys.exit(0)
