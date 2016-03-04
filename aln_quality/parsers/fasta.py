import logging
import os

from error_types import ParserError

_log = logging.getLogger(__name__)


def parse_fasta(aln_path, golden_ids):
    _log.info("Parsing FASTA: {}".format(aln_path))
    if not os.path.exists(aln_path):
        raise ParserError("FASTA file doesn't exist: {}".format(aln_path))
    with open(aln_path) as a:
        aln_file = a.read().splitlines()
    aln_dict = {}
    seq_id = ""
    for l in aln_file:
        if l.startswith(">"):
            if '|' in l:
                seq_id = l.lstrip('>').split('|')[0]
            else:
                seq_id = l.lstrip('>')
            if seq_id in aln_dict.keys():
                raise ParserError("Sequence {} is duplicated".format(seq_id))
            if seq_id in golden_ids:
                aln_dict[seq_id] = ""
        elif seq_id and seq_id in golden_ids:
            aln_dict[seq_id] += l.upper().rstrip("*")
    return aln_dict
