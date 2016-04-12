import os

from gold_standard_src.gold_standard.parsers.error_types import ParserError


def parse_3SSP(aln_path):
    if not os.path.exists(aln_path):
        raise ParserError("File not found: {}".format(aln_path))
    with open(aln_path) as a:
        infile = a.read().splitlines()
    seq_dict = {}
    for l in infile:
        seq_id = ''.join(l.split()[0:2])
        sequence = l.split()[2].replace('?', '-')
        seq_dict[seq_id] = sequence
    return seq_dict
