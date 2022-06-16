class ReadRecord:
    """Structure class that contains sequence reads information.
    srr_id: which SRR file it came from
    fragment_id: which fragment id was it
    """
    __slots__ = ("srr_id", "fragment_id", "sequence", "length")

    def __init__(self, srr_id, fragment_id, sequence, length):
        self.srr_id = srr_id
        self.fragment_id = fragment_id
        self.sequence = sequence
        self.length = length