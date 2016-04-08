import logging
import os

from gold_standard_src.gold_standard.num_seq import core_aln_to_num
from gold_standard_src.gold_standard.parsers.error_types import ParserError
from gold_standard_src.gold_standard.parsers.var_file import parse_var_file

_log = logging.getLogger(__name__)
fs = frozenset


def parse_gold_pairwise(gold_dir):
    _log.info("Getting gold standard alignments")
    if not os.path.exists(gold_dir):
        raise ParserError("No such directory: {}".format(gold_dir))
    # get all ".Var" files in the given directory
    var_list = [x for x in os.listdir(gold_dir) if x.endswith(".Var")]
    _log.info("Got %s var files", len(var_list))
    gold_alns = {}
    gold_ids = set()
    full_seq = {}
    for v in var_list:
        var = parse_var_file(os.path.join(gold_dir, v))
        gold_ids.update(var['ids'])
        aln_ids = fs(var['ids'])
        gold_alns[aln_ids] = var['alns']
        full_seq.update(var['full_seq'])
    _log.info("Finished parsing .Var files")
    return {
        'alns': gold_alns,
        'ids': gold_ids,
        'full_seq': full_seq
    }


def parse_gold_multi(gold_path, core_indexes):
    corvar = parse_var_file(gold_path)
    full_seq = get_full_seq(corvar)
    core_aln = remove_vars(corvar)
    num_aln = core_aln_to_num(core_aln, full_seq, core_indexes)
    return {'alns': num_aln, 'full_seq': full_seq,
            'ids': num_aln.keys()}


def remove_vars(corvar):
    # TODO: implement
    return corvar


def get_full_seq(corvar):
    # TODO: implement
    return corvar
