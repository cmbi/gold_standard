import os

from .error_types import ParserError


def parse_3SSP(aln_path):
    if not os.path.exists(aln_path):
        raise ParserError("File not found: {}".format(aln_path))
    with open(aln_path) as a:
        infile = a.read().splitlines()
    seq_dict = {}
    strcts_order = []
    for line in infile:
        if not line:
            # skip empty lines
            continue

        seq_id = ''.join(line.split()[0:2])
        seq_start = 2
        if len(seq_id) != 5:
            seq_id = line.split()[0]
            seq_start = 1
        strcts_order.append(seq_id)
        sequence = line.split()[seq_start:]
        seq_with_lower = [
            core[0].lower() + core[1:-1] + core[-1].lower() for core in sequence
        ]
        sequence = "".join(seq_with_lower)
        seq_dict[seq_id] = sequence

    return seq_dict, strcts_order
