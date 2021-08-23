import os
import sys
import pysam
from pybam.bamlib import BAM
from pybam.util import chunkstring

def zmw2fasta(bam_file, outdir):
    if not bam_file.endswith((".sam", ".bam")):
        sys.exit("Did you provide a SAM/BAM file?")

    state = 0
    alignments = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    for line in alignments:
        read = BAM(line)
        if state == 0:
            seq_lst = []
            current_zmw = read.zmw
            seq_lst.append((read.qname, read.qseq))
            state = 1 
        elif state == 1:
            if read.zmw == current_zmw:
                seq_lst.append((read.qname, read.qseq))
                state = 1
            else:
                state = 2

        if state == 2:
            # return
            o = open(os.path.join(outdir, "{}.fasta".format(current_zmw)), "w")
            for (qname, qseq) in seq_lst:
                o.write(">{}\n".format(qname))
                for chunk in chunkstring(qseq):
                    o.write("{}\n".format(chunk))
            o.close()
            # init
            seq_lst = []
            current_zmw = read.zmw
            seq_lst.append((read.qname, read.qseq))
            state = 1 
