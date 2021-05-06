#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"


# modules
import os
import sys
import logging
from pybam.fastq import bam2fastq
from pybam.parse_args import parse_args

def main():
    parser, options = parse_args(program_version=__version__)

    if options.sub == "bam2fastq":  # return first n lines of sequences
        bam2fastq(options.input, options.output)
    elif options.sub == "zmwlist":  # return first n lines of sequences
        bam2fastq(options.input, options.output)

    # elif options.sub == "filter": ## return first n lines of sequences
    #     hard_filter(options.input, options.number, options.output)
    # elif options.sub == "tricounts": ## returns counts based on trinucleotide context
    #     seq_statistics(options.input, options.ref, options.output)
    else:
        logging.warning("The subcommand does not exist!\n")
        parser.print_help()

if __name__ == "__main__":
    main()
    sys.exit(0)
