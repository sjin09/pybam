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
        required=True,
        type=argparse.FileType("w"),
        help=".fq, .fq.gz, .fastq or .fastq.gz file",
    )
    # subcommands: bam2fasta
    parser_head = subparsers.add_parser(
        "bam2fasta",
        help="converts BAM file to FASTA file",
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
        help=".fa, .f.gz, .fasta or .fasta.gz file",
    )
    # subcommands: zmwstat
    parser_head = subparsers.add_parser(
        "zmwstat",
        help="returns ZMW productivity statistics from CCS and subread BAM files",
    )
    parser_head.add_argument(
        "--ccs",
        type=str,
        required=True,
        help="PacBio circular consensus sequence (CCS) SAM/BAM or FASTQ (.fq, .fq.gz, fastq, fastq.gz) file",
    )
    parser_head.add_argument(
        "--subreads",
        type=str,
        required=True,
        help="PacBio subreads SAM/BAM file",
    )
    parser_head.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return the ZMW sequence statistics",
    )

    # subcommands: zmwlist
    parser_head = subparsers.add_parser(
        "zmwlist",
        help="returns zmwlist from BAM or FASTQ files",
    )
    parser_head.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="SAM/BAM or FASTQ file",
    )
    parser_head.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return list of zmws",
    )


    if len(arguments) == 0: # option length
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
