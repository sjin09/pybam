# modules
import sys
import argparse


# argparse
def parse_args(program_version, arguments=sys.argv[1:]):

    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="pybam parses SAM/BAM files and returns an output",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )

    # subcommands: bam2fastq
    subparsers = parser.add_subparsers(dest="sub")
    parser_head = subparsers.add_parser(
        "bam2fastq",
        help="converts BAM files to FASTQ files",
    )
    parser_head.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="SAM/BAM file",
    )
    parser_head.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help=".fq, .fq.gz, .fastq or .fastq.gz file",
    )

    # subcommands: subread_lengths
    parser_head = subparsers.add_parser(
        "subread_lengths",
        help="",
    )
    parser_head.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="SAM/BAM file",
    )
    parser_head.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return the ZMW and comma separated subread lengths",
    )

    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
