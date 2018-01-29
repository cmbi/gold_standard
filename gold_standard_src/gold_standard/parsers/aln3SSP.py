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
        if len(seq_id) != 5:
            seq_id = l.split()[0]
        strcts_order.append(seq_id)
        sequence = l.split()[2].replace('?', '-')
        seq_dict[seq_id] = sequence
    # remove gap-only sequences

    #for seq_id in seq_dict.keys():
        #if seq_dict[seq_id].count("-") == len(seq_dict[seq_id]):
            #del seq_dict[seq_id]
    return seq_dict, strcts_order
