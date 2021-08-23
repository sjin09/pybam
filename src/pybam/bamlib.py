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
        self.bq_int_lst = line.query_qualities
        self.bq_ascii = "".join([chr(_ + 33) for _ in line.query_qualities])
