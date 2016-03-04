import logging
import os

from error_types import ParserError
from var_file import parse_var_file

_log = logging.getLogger(__name__)
fs = frozenset


def parse_golden_alns(golden_dir):
    _log.info("Getting golden standard alignments")
    if not os.path.exists(golden_dir):
        raise ParserError("No such directory: {}".format(golden_dir))
    # get all ".Var" files in the given directory
    var_list = filter(lambda x: x.endswith(".Var"), os.listdir(golden_dir))
    _log.info("Got {} var files".format(len(var_list)))
    golden_alns = {}
    golden_ids = set()
    full_seq = {}
    for v in var_list:
        var = parse_var_file(os.path.join(golden_dir, v))
        golden_ids.update(var['ids'])
        aln_ids = fs(var['ids'])
        golden_alns[aln_ids] = var['aln']
        full_seq.update(var['full'])
    _log.info("Finished parsing .Var files")
    return golden_alns, golden_ids, full_seq
