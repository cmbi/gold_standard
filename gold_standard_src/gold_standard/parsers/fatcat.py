import logging
import os

from gold_standard_src.gold_standard.parsers.error_types import ParserError

_log = logging.getLogger(__name__)


def parse_fatcat(aln_path, golden_ids=None):
    _log.info("Parsing FATCAT output file: %s", aln_path)
    if not os.path.exists(aln_path):
        raise ParserError("File doesn't exist: {}".format(aln_path))
    with open(aln_path) as a:
        aln_file = a.read().splitlines()
    aln_dict = {}
    seq_id = ""
    for l in aln_file:
        seq_id = ''.join(l.split()[:2])
        if seq_id in golden_ids or not golden_ids:
            aln_dict[seq_id] = l.split()[2]
    return aln_dict
