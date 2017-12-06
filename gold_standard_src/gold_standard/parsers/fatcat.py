import logging
import os

from .parsers.error_types import ParserError

_log = logging.getLogger(__name__)


def parse_fatcat(aln_path, golden_ids=None):
    _log.info("Parsing FATCAT output file: %s", aln_path)
    if not os.path.exists(aln_path):
        raise ParserError("File doesn't exist: {}".format(aln_path))
    with open(aln_path) as a:
        aln_file = a.read().splitlines()
    aln_dict = {}
    seq_id = ""
    split_id = True
    for l in aln_file:
        seq_id = ''.join(l.split()[:2])
        if len(seq_id) != 5:
            seq_id = l.split()[0]
            split_id = False
        if seq_id in golden_ids or not golden_ids:
            if split_id:
                aln_dict[seq_id] = l.split()[2]
            else:
                aln_dict[seq_id] = l.split()[1]
    return aln_dict
