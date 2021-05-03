import pysam

class BAM:
    def __init__(self, line):
        # target
        self.tname = line.reference_name
        self.tstart = line.reference_start
        self.tend = line.reference_end
        # self.tcoord = "{}:{}-{}".format(self.tname, self.tstart, self.tend)
        # query
        self.qname = line.query_name
        self.qstart = line.query_alignment_start
        self.qend = line.query_alignment_end
        self.zmw = self.qname.split("/")[1]
        self.qseq = line.query_sequence
        self.qlen = len(line.query_sequence)
        self.bq_ascii = "".join([chr(_ + 33) for _ in line.query_qualities])
        # self.bq_int_lst = line.query_qualities
        # self.qcoord = "{}:{}-{}".format(self.zmw, self.qstart, self.qend)
        # alignment
        # self.mapq = str(line.mapping_quality)
        # self.orientation = "-" if line.is_reverse else "+"
        # self.alignment_length = line.query_alignment_length  # qend-qstart
        # self.alignment_fraction = str(self.alignment_length / self.qlen)
        # alignment tags
        # self.state = line.has_tag("cs")
        # self.cs_str = line.get_tag("cs") if line.has_tag("cs") else "."
        # self.alignment_type = line.get_tag("tp") if line.has_tag("tp") else "."
        # self.supplementary_alignments = (
        #     line.get_tag("SA:Z") if line.has_tag("SA:Z") else "."
        # )
