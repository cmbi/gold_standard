import os

from .error_types import ParserError


def parse_3SSP(aln_path):
    if not os.path.exists(aln_path):
        raise ParserError("File not found: {}".format(aln_path))
    with open(aln_path) as a:
        infile = a.read().splitlines()
    seq_dict = {}
    strcts_order = []
    for l in infile:
        seq_id = ''.join(l.split()[0:2])
        seq_start = 2
        if len(seq_id) != 5:
            seq_id = l.split()[0]
            seq_start = 1
        strcts_order.append(seq_id)
        # sequence = l.split()[2].replace('?', '-')
        sequence = l.split()[seq_start:]
        seq_with_lower = [
            core[0].lower() + core[1:-1] + core[-1].lower() for core in sequence
        ]
        sequence = "".join(seq_with_lower)
        seq_dict[seq_id] = sequence
    # remove gap-only sequences

    return seq_dict, strcts_order
